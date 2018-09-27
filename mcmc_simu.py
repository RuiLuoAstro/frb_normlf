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
bolupper=chkargv(sys.argv, '-upper')

#fname='simdat_ETG_43.0_-1.0.txt'
#fout='tsta'
frb_cat = lf.LoadSimuData(fname)

vFLUX = frb_cat['S']
vLOGFLUX = np.log10(vFLUX)
vDME=frb_cat['DMe']
#vDNU=frb_cat['dnu']
a1=time.clock()
#print dis.Normalize_constant_flux_Z(1,1000., 42., -1,37.)
a2=time.clock()
print a1, a2, a2-a1

vdnu = 1000.
logft = 0.
#mxlogls = np.loadtxt("./norm/simu/logls_simu.txt")
#myalpha = np.loadtxt("./norm/simu/alpha_simu.txt")
#mznorm = np.loadtxt("./norm/simu/norm_simu.txt")
#fintp = interpolate.interp2d(mxlogls[0,:], myalpha[:,0], mznorm, kind='linear')

def lik_fun(vpar):
    try:
        norm = dis.Norm1D(logft, vdnu, vpar[0], vpar[1], vpar[2])
        if norm==0:
            return -1e99
        res = np.sum(dis.log_distr_fdm(vdnu, vLOGFLUX, vDME, vpar[0], vpar[1], vpar[2], gtype=fgt)-np.log(norm))
        #print res
        return res
    except:
        print 'Numerical error: @', vpar
        return -1e99

#vdme=np.arange(0,3000,50)
#vflux=np.arange(0,4,0.1)
#mxdm, mxflux=np.meshgrid(vdme, vflux)
#vpar=[42., -1.0, 39.]
#mxres=np.zeros(mxdm.shape)

#n,m=mxdm.shape
#for i in range(n):
#    for j in range(m):
#        print i,j
#        dmv=mxdm[i,j]
#        fluxv=mxflux[i,j]
#        norm=dis.Normalize_1D(1,1000., vpar[0], vpar[1],vpar[2])
#        res=np.sum(dis.log_Distribution_flux_DM(1000, [fluxv], [dmv], vpar[0], 
#vpar[1], vpar[2], fgalaxy_type='ETG')-np.log(norm))
#        print res
#        mxres[i,j]=res

#np.savetxt('mxdm3', mxdm)
#np.savetxt('mxflux3', mxflux)
#np.savetxt('mxlik3', mxres)
#sys.exit(0)
if bolupper:
    vpara=np.array([40.0, -3.0, 1e35])
    vparb=np.array([46.0, 3.1, 1e40])
else:
    vpara=np.array([40.0, -3.0, 35.0])
    vparb=np.array([46.0, 3.1, 40.0])

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
                outputfiles_basename='mn_out/simu/'+fout)
