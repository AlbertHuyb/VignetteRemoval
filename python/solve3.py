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

# load pano file & images

import random
from merge import *
from arrayImage import *
# from numarray import *
import numpy as np
import sys
from functools import reduce

random.seed(2)

# pick random points in pano space
def rand2D(w,h):
    return (w*random.random(),h*random.random())

def nonzeroDRG(a):
    if len(a.shape)>1:
        print('nonzeroDRG only works for 1d arrays')
    nz = np.asarray(np.nonzero(a))
    if len(nz.shape)>len(a.shape):
        nz = np.reshape(nz,nz.shape[1:])
    return nz

def makeSamplePoints(panocam,panodict,npoints):
    # TODO use fromfunction
    p = np.transpose(np.array(list(map(lambda x:rand2D(*panocam.size),range(npoints)))))

    #print p.shape
    hitcounts = np.zeros(npoints)
    for image, camera in panodict.values():
        # project them into each image
        pc = panocam.cameraToCamera(camera,p)
        w,h = camera.size
        inside = np.logical_and.reduce([pc[0] >= 0,pc[0] < w,
                                     pc[1] >= 0,pc[1] < h])
        hitcounts += inside
    # save those which intersect at least 2 images
    #print 'before take:', p.shape
    p = np.take(p,nonzeroDRG(hitcounts >= 2),1)
    #p = compress(hitcounts >= 2,p,1)
    #print 'after take:', p.shape
    
    print(p.shape[1], "samples out of", npoints)
    # TODO reject samples with high gradients
    # TODO reject samples with saturated pixels
    return p

class VignData:
    error = 'VignData.error'
    
    def __init__(self,npixels,nimages):  # todo npixels is really npoints
        self.npixels = npixels
        self.nimages = nimages

        self.pindex = []
        self.iindex = []
        self.radius = []
        self.pixels = []
        
    def add(self,iindex,pindex,radius,pixels):
        if (len(pindex) != len(radius)):
            raise(error, "%d of radii (%d indices)" % (len(radius),len(pindex)))
        if (len(pindex) != len(pixels)):
            raise(error, "%d of pixels (%d indices)" % (len(pixels),len(pindex)))

        # TODO save iindex
        self.iindex += [iindex]*len(pindex)
        self.pindex += list(pindex)
        self.radius += list(radius)
        self.pixels += list(pixels)

    def complete(self):
        self.pindex = np.array(self.pindex)
        self.radius = np.array(self.radius)
        self.pixels = np.array(self.pixels)
        
    def __repr__(self):
        lens = [self.npixels] + map(len, [self.pindex,self.radius,self.pixels])
        return '<VignData: %d, %d %d %d>' % tuple(lens)

    def __len__(self):
        return len(self.pindex)

def makeVignData(panodict,panocam,p):
    data = VignData(p.shape[1],len(panodict))
    iindex = 0
    sortfiles = list(panodict.keys())
    sortfiles.sort()
    for filename in sortfiles:
        image,camera = panodict[filename]
        pc = np.asarray(panocam.cameraToCamera(camera,p))
        w,h = camera.size
        diag = np.sqrt(w*w+h*h)
        inside = np.logical_and.reduce([pc[0] >= 0,pc[0] < w,
                                     pc[1] >= 0,pc[1] < h])
        pindex = nonzeroDRG(inside)
        pc = np.take(pc,pindex,1)

        offset = pc - np.array([[w],[h]])/2
        radius = np.sqrt(sum(offset * offset))/(diag/2)
        #print radius.shape

        pixels = list(map(image.getpixel,list(map(tuple,np.transpose(pc)))))

        data.add(iindex,pindex,radius,pixels)
        iindex += 1

    data.complete()
    return data

class VignModel:
    error = 'VignModel.error'
    
    def __init__(self,curve=None,data=None,OneImageFlag=True):
        self.curve = curve
        self.dmax = None

        if data:
            self.data = data
            self.r2 = data.radius*data.radius
            self.r4 = self.r2*self.r2
            self.r6 = self.r4*self.r2

            if OneImageFlag:
                self.initCounts(OneImageFlag)
            else:
                self.initCounts()

            self.lamda = len(data)*10*10

            self.radiance = np.zeros((data.npixels,3),'d')
            
    def initCounts(self,OneImageFlag=False):
        data = self.data
        self.count = np.zeros(data.npixels)
        for i in range(len(data.pindex)):
            pindex = data.pindex[i]
            self.count[pindex] += 1
        # if minimum.reduce(self.count)<2:
        if (reduce(min,self.count)<2) and (OneImageFlag==False):
            raise(error, "some pixel only appears in 0 or 1 images")

    def getVignette(self,k):
        k2 = k[0]
        k4 = k[1]
        k6 = k[2]

        return np.ones(len(self.data)) + k2*self.r2 + k4*self.r4 + k6*self.r6

    def makeK(self,vk,rk,exposures):
        return vk + rk + exposures
    
    def splitK(self,k):
        k = np.asarray(k)
        vk = k[:3]
        rk = k[3:3+self.curve.size]
        exposures = k[3+self.curve.size:]
        return vk,rk,exposures
        
    def estRadiance(self,k):
        vk,rk,exposures = self.splitK(k)
        
        if len(exposures) != self.data.nimages:
            raise(self.error, "wrong # exposures (%d != %d)" % (len(exposures),self.data.nimages))
        v = self.getVignette(vk)
        self.curve.setCoeffs(rk,1)
        
        self.exposureat = np.take(exposures,self.data.iindex)
        estrad = self.curve.pixelsToLinear(self.data.pixels)
        # recover from vignette and exposure
        estrad /= (v*self.exposureat)[:,np.newaxis]

        data = self.data
        self.radiance = np.zeros((data.npixels,3),'d')
        #print self.data.pixels.shape, estrad.shape, self.radiance.shape
        for i in range(len(data.pindex)):
            pindex = data.pindex[i]
            self.radiance[pindex] += estrad[i]
        self.radiance /= self.count[:,np.newaxis]

    def estRadiance2(self,k):
        score = self.eval(k)
        print("previous score =", score)
        print("LINEAR RADIANCE ESTIMATION (FAST)")
        oldradiance = self.radiance.copy()
        self.estRadiance(k)       # initialize with linear
        nscore = self.eval(k)
        print("linear score =", nscore)
        if nscore < score:
            return                # good enough, do nothing

        self.radiance = oldradiance
        print("NONLINEAR RADIANCE ESTIMATION (SLOW)")
        from optimize import fmin
        for i in range(len(self.radiance)):
            if not i%500:
                print(i,"out of",len(self.radiance))
            for j in range(3):
                val = fmin(self.evalRad,(self.radiance[i,j],),(i,j),printmessg=0)
                self.radiance[i,j] = val[0]
        nnscore = self.eval(k)
        print("nonlinear score =", nnscore)
            
    def distance(self,a,b):
        diff = a-b
        if self.dmax:
            return sum(sum(np.minimum(diff*diff,self.dmax*self.dmax)))
        else:
            return sum(sum(diff*diff))
    
    def eval(self,k):
        vk,rk,exposures = self.splitK(k)
        if len(exposures) != self.data.nimages:
            raise(error, "wrong # exposures (%d != %d)" % (len(exposures),self.data.nimages))
        self.vign = self.getVignette(vk)
        # self.curve.setCoeffs(rk)
        self.curve.setCoeffs(rk,1)
        if hasattr(self.curve,'negLogLikelihood'):
            score = self.lamda * self.curve.negLogLikelihood(rk)
        else:
            score = 0
        
        radianceat = np.take(self.radiance,self.data.pindex,0)
        self.exposureat = np.take(exposures,self.data.iindex)
        recon = (self.vign*self.exposureat)[:,np.newaxis] * radianceat
        # recon = recon[:,np.newaxis]
        if self.curve:
            reconPixels = self.curve.linearToPixels(recon)
        else:
            reconPixels = recon
        diff = reconPixels - self.data.pixels
        score += sum(sum(diff*diff))

        #print 'eval', score
        return score

    def evalRad(self,rad,i,j):
        # same as eval, but with varying radiance instead of k
        # index gives the index of the point to be estimated
        # NOTE assumes that self.curve, self.vign, self.exposures already set

        # TODO find members of self.data.pindex that match index
        matches = np.nonzero(self.data.pindex == i)
        exposure = self.exposureat[matches]
        vign = self.vign[matches]
        pix = self.data.pixels[matches]
        diff = self.curve.linearToPixels(exposure*vign*rad) - pix[:,j]
        return sum(diff * diff)
        
    def optK(self,k):
        from optimize import fmin

        loopcounts = 0
        
        lastfval = 10e10
        fval = 10e9

        try:
            while fval < lastfval:
                lastfval = fval
                lastk = k
                self.estRadiance2(k)
                k,fval,warnflag = fmin(self.eval,k,fulloutput=1,printmessg=0)
                print(k,fval)
                loopcounts += 1
        except KeyboardInterrupt:
            pass

        print(loopcounts, "minimization loops")

        return lastk

    def makeVignArray(self,size,k):
        w,h = map(float,size)
        diag = np.sqrt(w*w+h*h)
        x = (np.array(range(int(w)))-w/2)/(diag/2)
        y = (np.array(range(int(h)))-h/2)/(diag/2)

        r2 = np.add.outer(x*x,y*y)
        r4 = r2*r2
        r6 = r2*r4

        o = np.ones(size)

        print(o.shape, r2.shape)
        
        return o + k[0]*r2 + k[1]*r4 + k[2]*r6

def drawSamples(p,infile,outfile):
    import Image
    im = Image.open(infile)

    import ImageDraw
    draw = ImageDraw.Draw(im)
    map(lambda x,d=draw:d.ellipse((x[0]-4,x[1]-4,
                                   x[0]+4,x[1]+4),fill=(255,0,0)),transpose(p))
    im.save(outfile)

if __name__ == "__main__":
    # imgdir = sys.argv[1]                  # 'small.ap'
    # npoints = int(sys.argv[2])            # 1000
    imgdir = 'images/senore/'
    npoints = 1000
    if len(sys.argv)>3:
        dmax = float(sys.argv[3])
    else:
        dmax = None

    kfile = imgdir + '/k.txt'
    emorfile = 'data/emor.txt'
    # pemorfile = 'data/pemor.txt'
    pemorfile = 'data/pemor.npy'

    panodict = parsePanoDir(imgdir)
    cameras = list(map(lambda x:x[1],panodict.values()))
    # Get the camera objects array

    from cylcam import makeCylCam
    panocam = makeCylCam(cameras)

    print("CHOOSING SAMPLES...")
    p = makeSamplePoints(panocam, panodict, npoints)

    unsolved = imgdir + '/unsolved.jpg'
    if (os.path.exists(unsolved)):
        drawSamples(p,unsolved,imgdir+'/samples.jpg')

    data = makeVignData(panodict,panocam,p)

    #from response import LookupCurve
    #curve = LookupCurve(crvfile)

    #from emor import EMoR
    print("LOADING RESPONSE MODEL...")
    from dorf import PEMoR
    curve = PEMoR(emorfile)
    curve.load(pemorfile)
    
    vm = VignModel(curve,data)

    if dmax:
        vm.dmax = dmax
        kinitfile = kfile
        kfile = imgdir + '/k.robust.txt'
        execfile(kinitfile)
        kinit = k
        del k
    else:
        kinit = vm.makeK([-.75,1,-.46],[0]*curve.size,[1]*len(cameras))

    print("KINIT = ", kinit)
    k = vm.optK(kinit)

    vk,rk,exposures = vm.splitK(k)

    print('# VK =', vk)
    print('# RK =', rk)
    print('# EXPOSURES =', exposures/exposures[0])
    print()

    fp = open(kfile,'w')
    fp.write('k = [ ')
    fp.write(', '.join(map(str,k)))
    fp.write(']')
    fp.close()
