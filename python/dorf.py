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

import string
# from numarray import *
import numpy as np
from emor import EMoR

dorffile = 'data/dorfCurves.txt'
emorfile = 'data/emor.txt'
pemorfile = 'data/pemor.txt'
pemorfile = 'data/pemor.npy'

def parseDorfFile(filename):
    #print 'reading file'
    fp = open(filename,'r')
    lines = fp.readlines()
    fp.close()

    #print 'parsing file'
    i = 0
    dorf = {}
    nlines = len(lines)
    while i < nlines:
        name = lines[i].strip();         i+= 1;
        curvetype = lines[i].strip();    i+= 1;
        I = []
        B = []
        i+=1
        while i < nlines and lines[i][0] in string.digits:
            I += map(float,lines[i].split())
            i+=1
        i+=1
        while i < nlines and lines[i][0] in string.digits:
            B += map(float,lines[i].split())
            i+=1
        
        dorf[name] = (type,np.asarray(I),np.asarray(B))

    return dorf

class PEMoR(EMoR):
    def setDorf(self,dorf):
        curves = map(lambda x: x[2], dorf.values())
        curves = np.array(list(curves))
        c = self.project(np.transpose(curves))
        cov = np.dot(c,np.transpose(c))

        from numpy.linalg import inv
        self.invsigma = inv(cov)

    def negLogLikelihood(self,c):
        ca = np.asarray(c)
        return np.dot(np.transpose(ca),np.dot(self.invsigma,ca))

    def save(self,filename):
        # fp = open(filename,'w')
        # fp.write(repr(self.invsigma))
        # fp.close()
        np.save(filename,self.invsigma)

    def load(self,filename):
        # fp = open(filename,'r')
        # s = "\n".join(fp.readlines())
        # fp.close()
        # self.invsigma = eval(s)
        self.invsigma = np.load(filename)
    
if __name__ == "__main__":
    dorf = parseDorfFile(dorffile)

    from emor import EMoR
    emor = EMoR(emorfile)

    pemor = PEMoR(emorfile)
    pemor.setDorf(dorf)
    pemor.save(pemorfile)

    pemor2 = PEMoR(emorfile)
    pemor2.load(pemorfile)
    nll = pemor2.negLogLikelihood
    
#    for curvetype, I, B in dorf.values():
#        c = emor.project(B)
#        print dot(transpose(c),dot(invsigma,c))

    print(nll([1,1,1,1,1]))
    print(nll([0.271685170126, 0.763604793495, -0.00019810334775,
               -0.0389643127791, 0.0720207625252]))
    print(nll([-1.78447991053, -2.42728052987, -0.487182472277,
               0.810791802423, 0.0929491177729]))
    print(nll([-1.28074878711, 0.236985072824, -0.170813711448,
               0.533237435281, -0.00836594147919]))
    
