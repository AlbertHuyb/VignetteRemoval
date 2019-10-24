import sys
import numpy as np
import cv2
from solve3 import rand2D, VignData, VignModel
from dorf import PEMoR

if __name__ == "__main__":
    # imgdir = sys.argv[1]                  # 'small.ap'
    # npoints = int(sys.argv[2])            # 1000
    imgdir = ['images/compound/camSnapshot_0.png','images/compound/camSnapshot_1.png']
    npoints = 5000

    init_kfile = 'images/compound/2000_k.txt'
    execfile(init_kfile)
    kfile = 'images/compound/5000_k.txt'
    emorfile = 'data/emor.txt'
    # pemorfile = 'data/pemor.txt'
    pemorfile = 'data/pemor.npy'
    homomatfile = 'data/1to0_1.5.npy'
    camera_num = 2

    OneImageFlag = True

    img = cv2.imread(imgdir[1])
    p = np.transpose(np.array(list(map(lambda x:rand2D(*(img.shape[0],img.shape[1])),range(npoints)))))
    p = p.astype(int)

    data = VignData(p.shape[1],camera_num)

    w = img.shape[0]
    h = img.shape[1]

    diag = np.sqrt(w*w+h*h)
    offset = p - np.array([[w],[h]])/2
    radius = np.sqrt(sum(offset * offset))/(diag/2)
    pindex = np.array(range(npoints))
    iindex = 0
    pixels = list(map(lambda x:img[x[0],x[1],:],list(map(tuple,np.transpose(p)))))

    data.add(iindex,pindex,radius,pixels)

    img2 = cv2.imread(imgdir[0])
    H = np.load(homomatfile)
    p2 = np.vstack((p,np.ones((1,p.shape[1]))))
    p2 = np.dot(H,np.array([p2[1],p2[0],p2[2]]))
    p2 = p2[:2,:]/p2[2]/1.5
    p2 = np.array([p2[1],p2[0]]).astype(int)

    offset2 = p2 - np.array([[w],[h]])/2
    radius2 = np.sqrt(sum(offset * offset))/(diag/2)
    pindex2 = np.array(range(npoints))
    iindex2 = 1
    pixels2 = list(map(lambda x:img2[x[0],x[1],:],list(map(tuple,np.transpose(p2)))))

    data.add(iindex2,pindex2,radius2,pixels2)
    data.complete()

    curve = PEMoR(emorfile)
    curve.load(pemorfile)
    
    vm = VignModel(curve,data)

    # kinit = vm.makeK([-.75,1,-.46],[0]*curve.size,[1]*camera_num)
    kinit = vm.makeK(k[:3],k[3:3+curve.size],k[3+curve.size:])
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
