#!/usr/bin/env python2.3
#########
#
# Copyright (c) 2005 Dan B Goldman
# 
# This file is part of the vignette-removal library.
#
# Vignette-removal is free software; you can redistribute it and/or modify
# it under the terms of the X11 Software License (see the LICENSE file
# for details).
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the X11
# License for more details.
#
#########

import os
import PIL.Image as Image
from functools import reduce

def panoExtent(basecam,cameras):
    # transform all camera's corners into frame of a
    w,h = basecam.size
    acorners = [(0,0),(0,h),(w,h),(w,0)]

    corners = map(lambda x,bc=basecam,cs=acorners: x.cameraToCamera(bc,cs),
                  cameras)
    corners = reduce(lambda a,b:a+b,corners)

    # find max extent of corners
    xs = map(lambda x:x[0],corners)
    ys = map(lambda x:x[1],corners)
    xmin = reduce(min,xs)
    xmax = reduce(max,xs)
    ymin = reduce(min,ys)
    ymax = reduce(max,ys)

    return (xmin, xmax, ymin, ymax)

def makePanoCam(basecam,cameras):
    '''a panoCam is the same as the basecam, but it has a different pixel size
    offset'''
    xmin, xmax, ymin, ymax = panoExtent(basecam,cameras)
    size = (int(xmax-xmin+0.5),int(ymax-ymin+0.5))
    offset = (-xmin,-ymin)
    from cameras import offsetCamera
    return offsetCamera(basecam,size,offset)

def makeMesh(panocam,othercam,ppblock=50):
    '''block by block transform'''
    nxblocks = panocam.size[0]/ppblock
    nyblocks = panocam.size[1]/ppblock
    xfracs = map(lambda x,nb=nxblocks: float(x)/nb,range(int(nxblocks) + 1))
    xfracs = list(xfracs)
    yfracs = map(lambda x,nb=nyblocks: float(x)/nb,range(int(nyblocks) + 1))
    yfracs = list(yfracs)
    mesh = []
    size = panocam.size
    for i in range(int(nxblocks)):
        bxspan = (int(size[0]*xfracs[i]),int(size[0]*xfracs[i+1]))
        for j in range(int(nyblocks)):
            byspan = (int(size[1]*yfracs[j]),int(size[1]*yfracs[j+1]))

            box = (bxspan[0],byspan[0],bxspan[1],byspan[1])
            boxcorners = [(box[0],box[1]),(box[0],box[3]),
                          (box[2],box[3]),(box[2],box[1])]
            pcorners = panocam.cameraToCamera(othercam,boxcorners)
            pcorners = reduce(lambda a,b:list(a)+list(b),pcorners)
            mesh.append((box,pcorners))
    return mesh

def parseSkipFile(file):
    skips = []
    fp = open(file)
    for line in fp.readlines():
        skips = skips + line.split()
    return map(int,skips)

def scaleCamera(image,camera,scale):
    ow, oh = image.size
    nw, nh = int(ow*scale), int(oh*scale)
    image = image.resize((nw,nh),Image.ANTIALIAS)
    camera.scale(nw,nh)
    return image,camera
    
def scaleCameras(panodict,scale):
    for filename in panodict.keys():
        image,camera = panodict[filename]
        panodict[filename] = scaleCamera(image,camera,scale)

def parsePanoDir(dir,scale=1):
    from cameras import parseCameraFile

    cameras = parseCameraFile(dir+'/pano.txt')

    dict = {}
    for camera in cameras:
        filename = dir + '/' + camera.filename.split('\\')[-1]
        #image = Image.open(filename).convert('RGBA')
        print('loading Image', filename)
        image = Image.open(filename)
        if scale != 1:
            print('...scaling by', scale)
            dict[filename] = scaleCamera(image,camera,scale)
        else:
            dict[filename] = (image,camera)

    skipfile = dir+'/skip.txt'
    if os.path.exists(skipfile):
        skips = parseSkipFile(skipfile)
        print("SKIPPING",skips)
        filenames = dict.keys()
        filenames.sort()
        for i in skips:
            del dict[filenames[i]]
        
    return dict

def makeTransformedImages(panodict,panocam):
    # reorder consistently:
    files = panodict.keys()
    files = list(files)
    files.sort()
    imgs = []
    for file in files:
        img,cam = panodict[file]
        mesh = makeMesh(panocam,cam)
        bnew = img.convert('RGBA').transform(panocam.size,Image.MESH,
                                             mesh,Image.BILINEAR)
        imgs.append(bnew)
    return imgs

def compositeAll(images):
    # create new image with max extent
    blend = Image.new('RGBA',images[0].size)
    for img in images:
        blend = Image.composite(img,blend,img)
    return blend

def makeBlendImage(panodict,panocam):
    # create new image with max extent
    blend = Image.new('RGBA',panocam.size)

    # reorder consistently:
    files = panodict.keys()
    files.sort()
    for file in files:
        img,cam = panodict[file]
        mesh = makeMesh(panocam,cam,30)
        bnew = img.transform(panocam.size,Image.MESH,mesh,Image.BILINEAR)
        blend = Image.composite(bnew,blend,bnew)

    return blend

if __name__ == '__main__':
    import sys, os
    #basecam = sys.argv[1]
    imgdir = sys.argv[1]
    outfile = sys.argv[2]
    if len(sys.argv) == 4:
        scale = float(sys.argv[3])
    else:
        scale = 1

    #imgdir = os.path.split(basecam)[0]

    panodict = parsePanoDir(imgdir,scale)

    #basecam = panodict[basecam][1]
    cameras = map(lambda x:x[1],panodict.values())

    from cylcam import makeCylCam
    panocam = makeCylCam(cameras)
    #panocam = makePanoCam(basecam,cameras)

    timgs = makeTransformedImages(panodict,panocam)

    outbase = '.'.join(outfile.split('.')[:-1])
    print(outbase)
    i = 0
    for img in timgs:
        filename = '.'.join([outbase,str(i),'png'])
        print(filename)
        img.save(filename)
        i += 1
    blend = compositeAll(timgs)
    blend.show()
    blend.save(outfile)

