#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pymultinest


def getargv(argv, key):
    for i in range(0,len(argv)):
        arg=argv[i]
        if (arg==key):
            return argv[i+1]

def chkargv(argv, key):
    for i in range(0,len(argv)):
        arg=argv[i]
        if (arg==key):
            return True

    return False

a = pymultinest.Analyzer(n_params = 8, outputfiles_basename=getargv(sys.argv, "-f"))
b = a.get_equal_weighted_posterior()
#allres=a.get_data()
#vlik=-allres[:,1]

mxlikli = b[:,-1]
mxchain = b[:,:-1]
#mxchain = allres
#rdat = np.array([[41, 45], [-3.0, 1.0], [35, 43]])
indmax = np.argmax(mxlikli)

alpha_best = mxchain[indmax,1]
logls_best = mxchain[indmax,0]

dat_best = np.array([alpha_best, logls_best])

rate = 1.0
indx = []
dat = mxchain
likli = mxlikli

row, col = dat.shape
if len(indx) == 0:
    indx = np.arange(0, col);

dat = dat.copy()
dat = dat[np.ix_(np.arange(int(row - row * rate), row), indx)]
likli = likli[np.ix_(np.arange(int(row - row * rate), row))]

par=np.array([])
parm=np.array([])
rangedat=np.array([])

if len(rangedat) == 0:
    rangedat = np.zeros((col, 2));
    for i in range(0, col):
        rangedat[i, 0] = np.min(dat[:, i]);
        rangedat[i, 1] = np.max(dat[:, i]);

npar = col
max_yticks = 4
levels = [0.68]

x = dat[:, 1]
y = dat[:, 0]

ind = (x < rangedat[1, 1]) & (x > rangedat[1, 0]) & (y < rangedat[0, 1]) & (y > rangedat[0, 0])
liklisub = likli[ind]

ngridx = 20
ngridy = 30

H, xedges, yedges = np.histogram2d(x, y, bins=(ngridx, ngridy), range=(rangedat[1, :], rangedat[0, :]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

xedges = (xedges[:-1] + xedges[1:]) / 2
yedges = (yedges[:-1] + yedges[1:]) / 2
mxx, mxy = np.meshgrid(xedges, yedges)

Ht = H.transpose()

indx, indy = np.meshgrid(np.arange(0, ngridx), np.arange(0, ngridy))
vmx2 = np.squeeze(np.reshape(mxx, (-1, 1)));
vmy2 = np.squeeze(np.reshape(mxy, (-1, 1)));
vm = np.squeeze(np.reshape(Ht, (-1, 1)));

vx = np.squeeze(np.reshape(indx, (-1, 1)));
vy = np.squeeze(np.reshape(indy, (-1, 1)));

vm2 = np.sort(vm)[::-1]
ix = np.argsort(vm, axis=0)[::-1]
ix = np.ix_(ix)
vx2 = vx[ix]
vy2 = vy[ix]
vmx2 = vmx2[ix]
vmy2 = vmy2[ix]

vm2 = np.cumsum(vm2 / np.sum(vm2))
cmxx2 = Ht
mxx2 = mxx
mxy2 = mxy

#print mxx2.shape
for ki in range(0, len(vm2)):
    mxx2[vy2[ki], vx2[ki]] = vmx2[ki]
    mxy2[vy2[ki], vx2[ki]] = vmy2[ki]
    cmxx2[vy2[ki], vx2[ki]] = vm2[ki]

conls = plt.contour(mxx2, mxy2, cmxx2, levels, colors='k');

condat = conls.collections[0].get_paths()[0]
vcondat = condat.vertices

if chkargv(sys.argv, '-out'):
	np.savetxt(getargv(sys.argv, "-ob"), dat_best)
	np.savetxt(getargv(sys.argv, "-oc"), vcondat)
else:
	print "No output, please try again."
