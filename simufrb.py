#!/bin/python

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
import time
import warnings
import sys

from mymodel import *

dis=AstroDistribution()
cos=Cosmology()
#frb_cat=LoadCatalogue()
#vBW=frb_cat['BW']
#vFLUX=frb_cat['S']
#vFLUX_thre=frb_cat['STD']
#vDM=frb_cat['DM']

#---------------Test code for Sampling N-dim----------------------------------
#def test_Gaussian(vpar):
#    v=np.exp(-0.5*(vpar[:,0]*vpar[:,0]+vpar[:,1]*vpar[:,1]))
#    return(v)
#
#res=SamplingND(zz, np.array([[-1,1],[-2,2]]), 1, 100000)
#np.savetxt('lala', res)
#exit(0)
#------------------------------------------------------------------------------

# Power-law index of FRB LF
alpha = float(getargv(sys.argv, '-alpha'))
# Cut-off luminosity
logls = float(getargv(sys.argv, '-logls'))
# Number of simulated FRBs
Ns=int(getargv(sys.argv, '-ns'))
#This is the threshold for selecting the source according to telescope
flux_thre=float(getargv(sys.argv, '-thre'))
# Bandwidth of the FRB survey
dnu=float(getargv(sys.argv, '-dnu'))
# galaxy type
fgt = getargv(sys.argv, '-type')

nt=0
res=np.zeros((Ns, 9))
Ns0=Ns
while (Ns>0):

    #Sampling the L
    vlogL=np.arange(38., 52., 12./10000)
    vlik=dis.Schechter_log(vlogL, 1, logls, alpha)
    vlogL=Sampling1D(vlogL, vlik,39, 51, Ns0)
    vlnEps=np.random.uniform(-np.log(2), 0, Ns0)
    vEps=np.exp(vlnEps)
    vZg=np.arange(0, 3.1, 3.1/10000)
    vlik=dis.Distribution_volume(vZg)
    vZ=Sampling1D(vZg, vlik, 0, 3.0, Ns0)
    vDMH0=np.arange(0, 5001., 5001./10000)
    vlik=dis.Distribution_HostGalaxyDM(vDMH0, fgalaxy_type=fgt)
    vDMH0=Sampling1D(vDMH0, vlik, 0, 5000, Ns0)
    #vDMH0=vDMH
    vDMH=vDMH0*np.sqrt(dis.SFR(vZ))/np.sqrt(dis.SFR(0))
    vDMS=np.random.uniform(0, 50, Ns0)
    vDMI=cos.DispersionMeasure_IGM(vZ)
    vDME=(vDMH+vDMS)/(1+vZ)+vDMI

    vFlux=vEps*cos.Luminosity_to_Flux(vZ, np.power(10., vlogL), dnu)
    nlen=len(vFlux[vFlux>flux_thre])
    if nlen>Ns:
        nlen=Ns

    res[nt:(nt+nlen),0]=vDME[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),1]=vFlux[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),2]=vlogL[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),3]=vZ[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),4]=vDMI[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),5]=vDMH[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),6]=vDMS[vFlux> flux_thre][0:nlen]
    res[nt:(nt+nlen),7]=dnu
    res[nt:(nt+nlen),8]=flux_thre
    nt=nt+nlen
    Ns=Ns-nlen
    print Ns, nlen

vDM=res[:,0]
vFlux=res[:,1]
vlogL=res[:,2]
vZ=res[:,3]
vDMe=res[:,4]
vDMh=res[:,5]
vDMs=res[:,6]
np.savetxt(getargv(sys.argv, '-o'), res, delimiter=' ',header="#DMe S logL Z DMi DMh DMs dnu thres", comments="")
