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

emorfile = 'data/emor.txt'

from solve3 import *
from merge import *
from functools import reduce
import operator

def removeVignettes(panodict,vm,k,outbase,curve=None):
    vk, rk, exposures = vm.splitK(k)

    vasize = (0,0)
    va = None

    print ('making response curve')
    curve.setCoeffs(rk,1)

    #exposures *= 10
    #exposures /= sum(exposures)/len(exposures)
    exposures /= pow(reduce(operator.mul,exposures),1.0/len(exposures))
    print (exposures)
    
    im = 0
    
    outdict = {}
    sortfiles = panodict.keys()

    sortfiles = list(sortfiles)
    sortfiles.sort()
#    for filename in sortfiles[:2]:
    for filename in sortfiles:
        print (filename)
        img, camera = panodict[filename]

        if vasize != camera.size:
            print ('making vignette array')
            vasize = camera.size
            va = vm.makeVignArray(vasize,vk)
            
        print ('image to array...')
        ia = image2array(img)
        
        if curve:
            ia = curve.pixelsToLinear(ia)
            print ('remove vignetting...')
            print (ia.shape, va.shape)
            devign = curve.linearToPixels(ia/(exposures[im]*va[...,np.newaxis]))
            # devign = curve.linearToPixels(ia/(va[...,np.newaxis]))
#            ia = curve.pixelsToLogR(ia)
#            print 'remove vignetting...'
#            devign = curve.logRToPixels(ia
#                                        - lva[...,NewAxis]
#                                        - log(exposures[im]))
    
        else:
            devign = ia / (exposures[im]*va[...,np.newaxis])

        print ('array to image...')
        dvi = array2image(devign).convert('RGBA')

        outfile = os.path.split(filename)[1]
        outfile = '.'.join([outbase] + outfile.split('.')[-2:])
        print ('saving', outfile)
        dvi.save(outfile,'png')

        outdict[outfile] = (dvi, camera)
        im += 1

    print ('done')
    return outdict

if __name__ == '__main__':
    import sys, os
    #basecam = sys.argv[1]
    # imgdir = sys.argv[1]
    # outfile = sys.argv[2]
    imgdir = 'images/senore/'
    outfile = 'test.png'
    outbase = '.'.join(outfile.split('.')[:-1])

    if len(sys.argv)>3:
        scale = float(sys.argv[3])
    else:
        scale = 1
        
    if len(sys.argv)>4:
        kfile = sys.argv[4]
    else:
        kfile = imgdir + '/k.old.txt'

    execfile(kfile)  # load k value
    
    panodict = parsePanoDir(imgdir,scale)
        
    #basecam = panodict[basecam][1]
    cameras = map(lambda x:x[1],panodict.values())
    cameras = list(cameras)

    from cylcam import makeCylCam
    panocam = makeCylCam(cameras)
    #panocam = makePanoCam(basecam,cameras)

    #vi = array2image(va)

    from emor import EMoR
    emor = EMoR(emorfile)
    vm = VignModel(emor)
    print (panodict.keys())
    panodict = removeVignettes(panodict,vm,k,outbase,emor)

    print (panodict.keys())
    
    timgs = makeTransformedImages(panodict,panocam)
    print (timgs)
    i = 0
    for img in timgs:
        filename = '.'.join([outbase,str(i),'png'])
        print (filename)
        img.save(filename)
        i += 1
    blend = compositeAll(timgs)
    # blend.show()
    blend.save(outfile)
