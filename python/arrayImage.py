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
UInt8 = 'UInt8'
import PIL.Image as Image

def array2image(a,scale=1):
    error = 'array2image.error'
    if len(a.shape) == 2:
        imtype = 'L'
    elif len(a.shape) == 3:
        if a.shape[2] == 3:
            imtype = 'RGB'
        elif a.shape[2] == 4:
            imtype = 'RGBA'
        else:
            raise (error, "wrong dimensions for image")
    else:
        raise (error, "wrong dimensions for image")

    data = np.minimum(a*scale,255.0).transpose(1,0,2)[:,:,::-1]
    # data = np.transpose(data.astype(UInt8),(1,0,2)).tostring()
    # size = a.shape[:2]
    # return Image.fromstring(imtype,size,data)
    return Image.fromarray(data.astype(np.uint8))

def image2array(i):
    # size = (i.size[1],i.size[0])
    # a = np.reshape(np.array(i.tostring(),UInt8),
    #             size+(len(i.mode),))
    a = np.asarray(i)
    a = a[:,:,::-1].transpose(1,0,2)
    return a
