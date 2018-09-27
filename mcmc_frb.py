#!/bin/python
import numpy as np
from scipy import interpolate
import time
import pymultinest
import warnings
from mymodel import *
import sys

dis = AstroDistribution()
cos = Cosmology()
lf = Loadfiles()

fname = getargv(sys.argv, '-f')
fout = getargv(sys.argv, '-o')
fgt = getargv(sys.argv, '-g')
bolremovehalo=chkargv(sys.argv, '-halo')
bolupper=chkargv(sys.argv, '-upper')
#fname='frb_cat.txt'
#fout='tst'
#fgt='ETG_NE2001'
specw = 1000.

frb_cat=lf.LoadCatalogue(fname)
vThre=lf.RMEq(frb_cat['SIG'], frb_cat['SEFD'], frb_cat['Npol'], frb_cat['BW'], frb_cat['W'])
vLOGFT=np.log10(vThre)
vFLUX = frb_cat['S']
vLOGFLUX = np.log10(vFLUX)


if fgt.find('NE2001')>=0:
    vDME=frb_cat['DME_NE2001']
else:
    vDME=frb_cat['DME_YMW16']

if bolremovehalo:
    vDME=vDME-30.

if fgt.find('ETG')>=0:
    fgt='ETG'

vDNU = np.ones(frb_cat['BW'].size)*specw
#print vBW



def lik_fun(vpar):
    try:
        norm = np.zeros(vLOGFT.shape)
        for i in range(0, len(norm)):
            logft = vLOGFT[i]
            norm[i] = dis.Norm1D(logft, specw, vpar[0], vpar[1], vpar[2])
        ind=norm<=0
        norm[ind]=1e-199
        res = np.sum(dis.log_distr_fdm(vDNU, vLOGFLUX, vDME, vpar[0], vpar[1], vpar[2], gtype=fgt)-np.log(norm))
        return res
    except:
        print 'Numerical error: @', vpar
        return -1e99

if bolupper:
    vpara=np.array([40.0, -3.0, 1e35])
    vparb=np.array([47.0, 3.1, 1e42])
else:
    vpara=np.array([40.0, -3.0, 35.0])
    vparb=np.array([47.0, 3.1, 42.0])

vpar_range=np.dstack((vpara.transpose(),vparb.transpose()))[0,:,:]

def myprior(cube, ndim, nparams):
    for i in range(ndim):
        cube[i] = vpar_range[i,0]+cube[i] *(vpar_range[i,1]-vpar_range[i,0])
    
    if bolupper:
        cube[2]=np.log10(cube[2])

def myloglike(cube, ndim, nparams):
    cube2=np.zeros(ndim)
    for i in range(0,ndim):
        cube2[i]=cube[i]

    res=lik_fun(cube2)
    return (res)

print '------------par range-----------'
print vpar_range
a1=time.clock()
if bolupper:
    vpar=vpara.copy()
    vpar[2]=np.log10(vpar[2])
    print myloglike(vpar, len(vpara),len(vpara))
else:
    print myloglike(vpara, len(vpara),len(vpara))

a2=time.clock()
print a1,a2
print "Start MCMC running ..."
# run MultiNest
pymultinest.run(myloglike, myprior, len(vpara),
                importance_nested_sampling = False,
                resume = False,
                verbose = True,
                sampling_efficiency = 'model',
                n_live_points = 1000,
                outputfiles_basename='mn_out/samp/'+fout)
