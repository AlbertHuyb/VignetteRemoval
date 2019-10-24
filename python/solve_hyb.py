import sys
import numpy as np
import cv2
from solve3 import rand2D, VignData, VignModel
from dorf import PEMoR

if __name__ == "__main__":
    # imgdir = sys.argv[1]                  # 'small.ap'
    # npoints = int(sys.argv[2])            # 1000
    imgdir = 'images/compound/camSnapshot_1.png'
    npoints = 5000

    emorfile = 'data/emor.txt'
    # pemorfile = 'data/pemor.txt'
    pemorfile = 'data/pemor.npy'

    OneImageFlag = True

    img = cv2.imread(imgdir)
    p = np.transpose(np.array(list(map(lambda x:rand2D(*(img.shape[0],img.shape[1])),range(npoints)))))
    p = p.astype(int)

    data = VignData(p.shape[1],1)

    w = img.shape[0]
    h = img.shape[1]

    diag = np.sqrt(w*w+h*h)
    offset = p - np.array([[w],[h]])/2
    radius = np.sqrt(sum(offset * offset))/(diag/2)
    pindex = np.array(range(npoints))
    iindex = np.zeros(pindex.shape)
    pixels = list(map(lambda x:img[x[0],x[1],:],list(map(tuple,np.transpose(p)))))

    data.add(iindex,pindex,radius,pixels)
    data.complete()

    curve = PEMoR(emorfile)
    curve.load(pemorfile)
    
    vm = VignModel(curve,data,OneImageFlag)

    kinit = vm.makeK([-.75,1,-.46],[0]*curve.size,[1]*1)

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
