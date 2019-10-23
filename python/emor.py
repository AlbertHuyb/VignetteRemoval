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

datafile = 'data/emor.txt'

# from numarray import *
import numpy as np

def searchsortedDRG(bin,values):
    return np.reshape(np.searchsorted(bin,np.ravel(values)),values.shape)

def parseFile(file):
    dict = {}
    array = None
    fp = open(file)
    for line in fp.readlines():
        if '=' in line:
            array = []
            name = line.split('=')[0].strip()
            dict[name] = array
        else:
            for word in line.split():
                array.append(float(word))
    fp.close()
    return dict

def lerp(xm,xM,x,ym,yM):
    return np.where(xM == xm,
                 ym,
                 (x-xm)/(xM-xm)*(yM-ym) + ym)

def makeMonotonic(a):
    last = a[0]
    for i in range(1,len(a)):
        a[i] = max(last,a[i])
        last = a[i]
    return a

def isMonotonic(a):
    return logical_and.reduce(a[1:] - a[:-1] >= 0)
        
class EMoR:
    def __init__(self,file,trunc=5):
        d = parseFile(file)
        self.x = np.asarray(d['E'])
        self.mean = np.asarray(d['f0'])
        self.basis = []
        self.size = trunc
        for i in range(trunc):
            self.basis.append(np.asarray(d['h(%d)' % (i+1)]))
        self.basis = np.transpose(self.basis)
        self.setCoeffs([0]*self.size)

    def project(self,curves):
        from numpy.linalg import lstsq as lls
        zmc = curves-self.mean[:,np.newaxis]
        x, ss, rank, sv = lls(self.basis,zmc)
        return x
        
    def setCoeffs(self,v,makeMono=0):
        self.rc = self.mean + np.dot(self.basis,v)
        if makeMono:
            self.rc = makeMonotonic(self.rc)
        return self.rc

    def isMonotonic(self):
        return isMonotonic(self.rc)

    def eval(self,x):
        # linear interpolation
        i = np.minimum(searchsortedDRG(self.x,np.array(x)),len(self.x)-2)
        return np.clip(
            lerp(self.x[i],self.x[i+1],x,
                 self.rc[i],self.rc[i+1]),
            0,1)

    def pixelsToLinear(self,a):
        fp = a/255.
        i = np.minimum(searchsortedDRG(self.rc,fp),len(self.rc)-2)
        return np.clip(
            lerp(self.rc[i],self.rc[i+1],fp,
                 self.x[i],self.x[i+1]),
                 0,1)
        
    def linearToPixels(self,a):
        return 255*self.eval(a)
        

if __name__ == '__main__':
    import sys
    coeffs = map(float,sys.argv[1:])
    emor = EMoR(datafile,len(coeffs))
    # emor = EMoR(datafile)
    emor.setCoeffs(coeffs,1)

    if 0:
        l = [-.1,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1]
        #print emor.eval(l)
        p = emor.linearToPixels(l)
        #print p
        l = emor.pixelsToLinear(p)
        #print l
    else:
        print(emor.pixelsToLinear(array([0])))
        for i in range(len(emor.x)):
            print(emor.x[i],emor.rc[i])
        
