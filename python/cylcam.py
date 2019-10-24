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

from types import *
from math import *
# from numarray import *
# from numarray.linear_algebra import *
import numpy as np
import numpy.linalg as linalg
from functools import reduce

class CylCam:
    def __init__(self, size, thetarange, yrange):
        self.size = size
        offset = np.array([[1,0,-thetarange[0]],
                        [0,1,-yrange[0]],
                        [0,0,1]])
        scale = np.array([[size[0]/(thetarange[1]-thetarange[0]),0,0],
                       [0,size[1]/(yrange[1]-yrange[0]),0],
                       [0,0,1]])

        self.proj = np.dot(scale,offset)
        self.invproj = linalg.inv(self.proj)

    def worldToCamera(self,p):
        theta = np.arctan2(p[1],p[2])
        y = p[0]/np.sqrt(p[1]*p[1]+p[2]*p[2])
        c = np.dot(self.proj,[theta,y,1])
        return np.array([c[0],c[1]])

    def cameraToWorld(self,p):
        if len(p) == 2:
            tp0 = type(p[0])
            if isinstance(p[0],(float,int)):
                p = np.array([p[0],p[1],1.0])
            elif isinstance(p[0],np.ndarray):
                #print p.shape, p[0].shape, p[1].shape
                p = np.array([p[0],p[1],np.ones(p[0].shape)])
        else:
            p = np.array(p)
        ty = np.dot(self.invproj,p)
        theta = ty[0]/ty[2]
        y = ty[1]/ty[2]
        return np.array([y,np.sin(theta),np.cos(theta)])

    def cameraToCamera(self,ocam,p):
        if isinstance(p,list):
            return list(map(lambda x,ocam=ocam:self.cameraToCamera(ocam,x),p))
        else:
            return ocam.worldToCamera(self.cameraToWorld(p))

def yExtent(cyl,p0,p1):
    p0s = np.array([p0[1],p0[0],p0[2]])
    p1s = np.array([p1[1],p1[0],p1[2]])
    d = p1s-p0s
    denom = -d[0]*d[1]*p0s[0]+d[0]*d[0]*p0s[1]+d[2]*(d[2]*p0s[1]-d[1]*p0s[2])
    if denom == 0:
        pc = p0
    else:
        num = -p0s[1]*(d[0]*p0s[0]+d[2]*p0s[2])+d[1]*(p0s[0]*p0s[0]*p0s[1]*p0s[1])
        t = num/denom
        t = min(max(0,t),1)
        #print 't',t
        pc = p0+(p1-p0)*t
    #pc = (p1+p0)/2
    pts = map(lambda x,c=cyl:c.worldToCamera(x),[p0,pc,p1])

    #print pts
    ys = map(lambda x:x[1],pts)
    ys = list(ys)
    ymin = reduce(min,ys)
    ymax = reduce(max,ys)
    #print ymin, ymax
    return ymin,ymax

def cylExtent(cameras):
    cyl = CylCam((2*pi,1),(0,2*pi),(0,1))
    w,h = cameras[0].size
    corners = [(0,0),(0,h),(w,h),(w,0)]
    ccorners = map(lambda x,cyl=cyl,cnr=corners: x.cameraToCamera(cyl,cnr),
                  cameras)
    ccorners_result = list(ccorners)
    ccorners = reduce(lambda a,b:a+b,ccorners_result)

    ts = map(lambda x:x[0],ccorners)
    ts = list(ts)
    tmin = reduce(min,ts)
    tmax = reduce(max,ts)

    ymin = 10e6
    ymax = -10e6
    for camera in cameras:
        wcorners = map(lambda x,c=camera:c.cameraToWorld(x),corners)
        wcorners = list(wcorners)
        for i in range(len(wcorners)):
            p0 = np.array(wcorners[i])
            p1 = np.array(wcorners[(i+1)%4])
            #print 'edge',i
            ymini,ymaxi = yExtent(cyl,p0,p1)
            ymin = min(ymini,ymin)
            ymax = max(ymaxi,ymax)

    return (tmin, tmax, ymin, ymax)

def makeCylCam(cameras,scale=1):
    tmin, tmax, ymin, ymax = cylExtent(cameras)
    print (tmin, tmax, ymin, ymax)
    scale *= cameras[0].focal
#    return CylCam((800,300),(-pi/4,pi/4),(0,1))
    size = (int(scale*(tmax-tmin)),int(scale*(ymax-ymin)))
    print ('making cylcam of size', size)
    return CylCam(size,(tmin,tmax),(ymin,ymax))
