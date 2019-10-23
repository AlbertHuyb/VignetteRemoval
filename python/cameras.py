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

# from numarray import *
import numpy as np
import numpy.linalg as linalg
from types import *

class Camera:
    def __init__(self, filename, size, T, R, focal):
        self.filename = filename.strip()
        self.size = size
        self.T = np.array(T)
        self.R = np.array(R)
        self.focal = focal

        self.valid = np.logical_or.reduce(np.ravel(self.R != 0))
        if not self.valid:
            return

        # Projection = T * K(f) * R
        self.K = np.array([[focal, 0, 0],[0, focal, 0], [0, 0, 1]])

        swapxy = np.array([[0,1,0],[1,0,0],[0,0,1]])
        
        self.proj = linalg.multi_dot([swapxy,self.T,self.K,self.R])
        self.invproj = linalg.inv(self.proj)

    def scale(self, nw, nh):
        sx = float(nw)/self.size[0]
        sy = float(nh)/self.size[1]
        self.size = (nw,nh)
        self.focal *= (sx+sy)/2

        self.proj = np.dot(np.array([[sx,0,0],[0,sy,0],[0,0,1]]),self.proj)
        self.invproj = linalg.inv(self.proj)
        
    def __repr__(self):
        return "<Camera \"%s\": %s>" % (self.filename, self.proj)

    def worldToCamera(self,p):
        c = np.dot(self.proj,p)
        if (len(c.shape) == 1):
            if c[2]>0:
                return np.array([c[0]/c[2],c[1]/c[2]])
            else:
                return np.array([1e6,1e6])
        else:
            selector = np.repeat((c[2]>0)[np.newaxis,:],2,0)
            n = c.shape[1]
            default = np.repeat(np.array([1e6,1e6])[:,np.newaxis],n,1)
            return np.where(selector,np.array([c[0]/c[2],c[1]/c[2]]),default)
#        if c[2]>0:
#            return array([c[0]/c[2],c[1]/c[2]])
#        else:
#            return array([1e6,1e6])

    def cameraToWorld(self,p):
        #print 'cw'
        if len(p) == 2:
            tp0 = type(p[0])
            if isinstance(p[0],(float,int)):
                p = np.array([p[0],p[1],1.0])
            elif isinstance(p[0],np.ndarray):
                #print p.shape, p[0].shape, p[1].shape
                p = np.array([p[0],p[1],np.ones(p[0].shape)])
        else:
            p = np.array(p)
        #print p.shape
        return np.dot(self.invproj,p)

    def cameraToCamera(self,ocam,p):
        # if type(p) == ListType:
        if isinstance(p,list):
            return list(map(lambda x,ocam=ocam:self.cameraToCamera(ocam,x),p))
        else:
            return ocam.worldToCamera(self.cameraToWorld(p))

def offsetCamera(ocam,size,offset):
    cam = Camera('',size,ocam.T,ocam.R,ocam.focal)

    trans = array([[1,0,offset[0]],[0,1,offset[1]],[0,0,1]])
    cam.proj = dot(trans,cam.proj)
    cam.invproj = inverse(cam.proj)
    return cam
    
def parseCamera(fp):
    filename = fp.readline()
    if not filename:
        raise EOFError
    
    size = map(int,fp.readline().split())
    size = list(size)
    fp.readline()

    T = []
    for i in range(3):
        T.append(list(map(float,fp.readline().split())))
    fp.readline()

    R = []
    for i in range(3):
        R.append(list(map(float,fp.readline().split())))
    fp.readline()

    focal = float(fp.readline())
    fp.readline()

    retval = Camera(filename, size, T, R, focal)

    if retval.valid:
        return retval
    else:
        return None
    
def parseCameraFile(file):
    '''return a list of transforms'''
    fp = open(file)
    records = []
    while 1:
        try:
            cam = parseCamera(fp)
            if cam:
                records.append(cam)
        except EOFError:
            break
        
    return records

if __name__ == '__main__':
    panofile = 'small/pano.txt'
    pixel = [161,20]
    #panofile = '2005_01_11/manual.1/pano.txt'
    #pixel = [956, 902]
    #pixel = [643, 82]

    cameras = parseCameraFile(panofile)

    c1 = cameras[0]
    c2 = cameras[2]

    newpixel = c2.worldToCamera(c1.cameraToWorld(pixel))
    print (newpixel)
