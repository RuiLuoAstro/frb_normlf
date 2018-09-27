#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
from mpl_toolkits.mplot3d import Axes3D
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

def plot2dposterier_withconf(dat, likli, indx=[], labels=[], rate=0.5, 
        par=np.array([]),
                    parm=np.array([]), rangedat=np.array([]), levels=[0.68]):
    '''Plot a 2D posterier with likelihood burning curve
    dat: a n x m numpy matrix, each data point is a row
    indx: a array of integer indicating wich column (parameter) to plot
    labels: the LaTex label of the given parameter
    rate: the percentage of data to plot, e.g. 0.3 means plot data from 70% to
    the end.
    '''

    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font', **{'family': 'serif', 'serif': ['Palatino'], 'size': 8})
    #rc('text', usetex=True)
    row, col = dat.shape;
    if len(indx) == 0:
        indx = np.arange(0, col);

    if len(labels) == 0:
        labels = ['a'];
        labels = labels * (col);
    dat=dat.copy()
    dat = dat[np.ix_(np.arange(int(row - row * rate), row), indx)];
    likli = likli[np.ix_(np.arange(int(row - row * rate), row))];
    row, col = dat.shape;
    #max_yticks = 3
    max_yticks = 4
    if len(rangedat) == 0:
        rangedat = np.zeros((col, 2));
        for i in range(0, col):
            rangedat[i, 0] = np.min(dat[:, i]);
            rangedat[i, 1] = np.max(dat[:, i]);

    npar = col;

    #index of y
    for vari in range(0, npar):
        # plot the histogram
        #ax = plt.subplot2grid((npar, npar), (vari, vari))

        ax=plt.subplot(npar, npar, vari*npar+vari+1)
        ind = ((dat[:, vari] < rangedat[vari, 1]) & (dat[:, vari] > rangedat[vari, 0]))
        n, bins, patches = plt.hist(dat[ind, vari], 100, normed=True, \
                                    histtype='stepfilled', range=(rangedat[vari, 0], rangedat[vari, 1]))
        plt.setp(patches, 'facecolor', 'y', 'alpha', 0.75)
        yloc = plt.MaxNLocator(max_yticks)
        ax.xaxis.set_major_locator(yloc)
        xloc = plt.MaxNLocator(max_yticks)
        ax.yaxis.set_major_locator(xloc)
        plt.xlabel(labels[vari])
        if len(levels)>0:
            
            for ls in [0.68, 0.95]:
                histc, bin_edges = np.histogram(dat[:, vari], density=True, bins=100)
                vom=(bin_edges[:-1]+bin_edges[1:])*0.5
                ind=np.argsort(histc)
                ind=ind[::-1]
                v=histc.copy()
                v[0]=histc[ind[0]]
                for i in range(1, len(ind)):
                    v[i] = v[i-1]+histc[ind[i]]

                v=v/float(np.sum(histc))
                thre= histc[ind[ (v[:-1] < ls) & (v[1:]>=ls)][0]]

                lv=np.min(vom[histc>=thre])
                rv=np.max(vom[histc>=thre])
                print vari,'-th parameter', 'sigma=', ls, 'lv=', lv, 'rv=', rv, thre
                
                #if vari == 2 and ls == 0.95:
                #    plt.plot([rv, rv], [0, max(n)], ls='dashed', color='k', linewidth=1)
                if vari != 2:
                    plt.plot([lv, lv], [0, max(n)], ls='dashed', color='k', linewidth=1)
                    plt.plot([rv, rv], [0, max(n)], ls='dashed', color='k', linewidth=1)
                    print rv-par[vari], lv-par[vari]
                elif ls == 0.95:
                    print rv
                    plt.plot([rv, rv], [0, max(n)], ls='solid', color='k', linewidth=2)
            
        if par.size > 0 and vari != 2:
            plt.plot([par[vari], par[vari]], [0, max(n)], ls='solid',
                     color='k', linewidth=2)
        
        if len(parm) > 0:
            for i in range(len(parm)):
                parm0=parm[i]
                plt.plot([parm0[vari], parm0[vari]], [0, max(n)], ls='dashed',
                     color='k', linewidth=2)
        
        plt.xlim(rangedat[vari, :])
        #index of x
        for varj in range(vari + 1, npar):

            #ax = plt.subplot2grid((npar, npar), (vari, varj))
            ax=plt.subplot(npar, npar, (vari)*npar+varj+1)
            x = dat[:, varj]
            y = dat[:, vari]

            ind = (
            (x < rangedat[varj, 1]) & (x > rangedat[varj, 0]) & (y < rangedat[vari, 1]) & (y > rangedat[vari, 0]))
            liklisub = likli[ind]
            x = dat[ind, varj]
            y = dat[ind, vari]
            ngridx = 20
            ngridy = 30
            #generate 2D histogram
            H, xedges, yedges = np.histogram2d(x, y, bins=(ngridx, ngridy),
                                               range=(rangedat[varj, :], rangedat[vari, :]))
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            #hisogram is row leading need transpose to plot with contourf
            H = H.transpose()
            #get the center of bin
            xedges = (xedges[:-1] + xedges[1:]) / 2
            yedges = (yedges[:-1] + yedges[1:]) / 2
            mxx, mxy = np.meshgrid(xedges, yedges)
            plt.contourf(mxx, mxy, H, 100, cmap=plt.cm.summer_r);
            plt.title("%(aa)s-%(bb)s" % {'aa': labels[vari], 'bb': labels[varj]})
            #Prepare to count the confidence level
            indx, indy = np.meshgrid(np.arange(0, ngridx), np.arange(0, ngridy))
            vmx2 = np.squeeze(np.reshape(mxx, (-1, 1)));
            vmy2 = np.squeeze(np.reshape(mxy, (-1, 1)));
            vm = np.squeeze(np.reshape(H, (-1, 1)));

            vx = np.squeeze(np.reshape(indx, (-1, 1)));
            vy = np.squeeze(np.reshape(indy, (-1, 1)));
            #Sort according to likelihood
            vm2 = np.sort(vm)[::-1];
            #Get index and convert everything to array
            ix = np.argsort(vm, axis=0)[::-1];
            ix = np.ix_(ix);
            vx2 = vx[ix];
            vy2 = vy[ix];
            vmx2 = vmx2[ix];
            vmy2 = vmy2[ix];
            #get the cumulative of the hisogram
            vm2 = np.cumsum(vm2 / np.sum(vm2));
            cmxx2 = H;
            mxx2 = mxx;
            mxy2 = mxy;
            #Form 2D cumulative plot
            for ki in range(0, len(vm2)):
                mxx2[vy2[ki], vx2[ki]] = vmx2[ki];
                mxy2[vy2[ki], vx2[ki]] = vmy2[ki];
                cmxx2[vy2[ki], vx2[ki]] = vm2[ki];
            conls = plt.contour(mxx2, mxy2, cmxx2, levels, colors='k');
            plt.clabel(conls, inline=1, fontsize=10)
            ax.get_yaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)
            yloc = plt.MaxNLocator(max_yticks)
            ax.xaxis.set_major_locator(yloc)
            xloc = plt.MaxNLocator(max_yticks)
            ax.yaxis.set_major_locator(xloc)

a = pymultinest.Analyzer(n_params = 8, outputfiles_basename=getargv(sys.argv, "-f"))
b = a.get_equal_weighted_posterior()
#allres=a.get_data()
#vlik=-allres[:,1]

vlik=b[:,-1]
allres=b[:,:-1]
mxchain = allres
rdat = np.array([[41, 46], [-3.0, 3.0], [38, 43]])
vpar = mxchain[np.argmax(vlik),:]
#print vlik
print vpar
mpl.rcParams['font.size']=18
mpl.rcParams['xtick.labelsize']=16
mpl.rcParams['ytick.labelsize']=16
mpl.rcParams['axes.labelsize']=20


if chkargv(sys.argv, '-o'):
	plt.figure(figsize=(10,7))
	plot2dposterier_withconf(mxchain, vlik, par=vpar, indx=[], rangedat=rdat, rate=1.0, levels=[0.68, 0.95], labels=[r'$\log L^*$', r'$\alpha$', r'$\log L_0$'])
	plt.suptitle(getargv(sys.argv, "-title"))
	plt.savefig(getargv(sys.argv, "-o"))

else:
	plot2dposterier_withconf(mxchain, vlik, par=vpar, indx=[], rangedat=rdat, rate=1.0, levels=[0.68, 0.95], labels=[r'$\log L^*$', r'$\alpha$', r'$\log L_0$'])
	plt.suptitle(getargv(sys.argv, "-title"))
	plt.show()
