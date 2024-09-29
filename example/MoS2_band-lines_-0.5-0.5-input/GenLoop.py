#!/usr/bin/env python
import os
import numpy as np
import shutil

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
n=40
for i in range(n):  #the mesh for k line of that wcc are calculated
    filename = 'file%d'%(i+1)
    workdir = os.path.join(MODULE_DIR, filename)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    os.chdir(workdir)

    shutil.copy(os.path.join(MODULE_DIR,'INCAR.band'),workdir)
    os.rename('INCAR.band','INCAR')
#    shutil.copy(os.path.join(MODULE_DIR,'INCAR.song.band'),workdir)
#    os.rename('INCAR.song.band','INCAR.song')
    shutil.copy(os.path.join(MODULE_DIR,'POSCAR'),workdir)
    shutil.copy(os.path.join(MODULE_DIR,'POTCAR'),workdir)
#    shutil.copy(os.path.join(MODULE_DIR,'CHGCAR'),workdir)
#    shutil.copy(os.path.join(MODULE_DIR,'CHG'),workdir)
    #
    os.system('ln -s '+os.path.join(MODULE_DIR,'CHGCAR')+' CHGCAR')
    os.system('ln -s '+os.path.join(MODULE_DIR,'CHG')+' CHG')
    
    kfiles = os.path.join(workdir,'KPOINTS')
    f = open(kfiles, 'w')
    f.write('K-points\n')
    f.write('%d\n'%n)
    f.write('rec\n')
    kx = -0.5+i*(1.0/n)   
    for j in range(n):
        ky = -0.5+j*(1.0/n)
        kz = 0
        f.write('%12.6f  %12.6f   %12.6f    1\n' %(kx, ky, kz))
    f.close()



#    do ik2=1, Nk2 ! the mesh for k line of that wcc are calculated
#        theta2= (ik2- 1d0)/(Nk2- 1d0)* pi
#        if (ik2== 1) theta2= (ik2- 1d0+ 0.10)/(Nk2- 1d0)* pi  ! avoid the North pole
#        if (ik2== Nk2) theta2= (ik2- 1d0- 0.10)/(Nk2- 1d0)* pi  ! avoid the south pole
#        do ik1=1, Nk1  ! the mesh along the integration direction
#           theta1= (ik1- 1d0)/Nk1* 2d0* pi
#           r_para= r0* sin(theta2)
#           kpoints(1, ik1, ik2)= k0(1)+ r_para* cos (theta1)
#           kpoints(2, ik1, ik2)= k0(2)+ r_para* sin (theta1)
#           kpoints(3, ik1, ik2)= k0(3)+ kr0* cos (theta2)
#        enddo
#     enddo
#irot  :   4   #m001
# --------------------------------------------------------------------
# isymop:   1   0   0
#           0   1   0
#           0   0  -1
# 
# gtrans:     0.0000000     0.0000000     0.0000000
# 
# ptrans:     0.0000000     0.0000000     0.0000000
# rotmap:
# (   1->   1)  (   2->   2)  (   3->   3)  (   4->   4)  (   5->   5) 
# (   6->   6)  (   7->   7)  (   8->   8)  (   9->   9) 
#
