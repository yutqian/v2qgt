#!/usr/bin/env python
import os
import numpy as np
import shutil

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
n=101
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
    #os.system('ln -s '+os.path.join(MODULE_DIR,'CHGCAR')+' CHGCAR')
    #os.system('ln -s '+os.path.join(MODULE_DIR,'CHG')+' CHG')
    
    kfiles = os.path.join(workdir,'KPOINTS.tmp')
    f = open(kfiles, 'w')
    #f.write('K-points\n')
    #f.write('%d\n'%n)
    #f.write('rec\n')

    kx = -0.17+0.34*i*(1.0/n)
    for j in range(n):
        ky =           -0.13 +0.26*j*(1.0/n)
        kz =  0.5+0.5*(-0.13 +0.26*j*(1.0/n))
        f.write('%12.6f  %12.6f   %12.6f    0\n' %(kx, ky, kz))
    f.close()

