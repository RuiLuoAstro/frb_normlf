import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import scipy.special as spf

class Cosmology:
    def __init__(self, omegam=0.308, omegal=0.692, omegab=0.0484):
        self.Omega_m = omegam
        self.Omega_b = omegab
        self.Omega_L = omegal
        self.c = 2.9979245800e10
        self.pc2cm = 3.08567758e18
        self.km2cm = 1e5
        self.Mpc2cm = 3.08567758e24
        self.Jy2CGS = 1e-23
        self.MHz2Hz = 1e6
        self.Jyms2CGS = 1e-26
        self.h0 = 0.6781
        self.H0 = self.h0 * 100 * self.km2cm / self.Mpc2cm
        self.Rhoc = 1.88 * self.h0 * self.h0 * 1e-29    #The critical density of universe in gram/cm^3
        self.Nc = self.Rhoc / 1.6726e-24       #The number density of universe in Hydrogen atom, in units of 1/cm^3
        self.f_IGM = 0.83
        self.z0 = 0.8

        self.vz = np.arange(-7, 6, 0.03)
        self.vz = np.power(10., self.vz)
        self.vz[0] = 0
        func = lambda z: self.c / self.H(z)
        self.vd = [ integrate.quad(func, 0, zv)[0] for zv in self.vz ]
        self.cd_interp = interp1d(self.vz, self.vd)
                
        func2 = lambda z: (1+z) * self.c / self.H(z) * self.Omega_b * self.Nc
        self.vdm = [ integrate.quad(func2, 0, zv)[0] for zv in self.vz ]
        self.dmigm = interp1d(self.vz, self.vdm)
        
        self.vld = self.Luminosity_Distance(self.vz)
        self.Ld2z = interp1d(self.vld, self.vz)
        self.Ld1 = self.Luminosity_Distance(1.0)
        self.Cd1 = self.Comoving_Distance(1.0)
        
    def E(self, z):
        """
        logarithmic time derivative of scale factor
        """
        return np.sqrt(self.Omega_m * np.power(1 + z, 3.0) + self.Omega_L)
    
    def H(self, z):
        """
        calculate the Hubble ratio
        z: redshifts
        """
        return self.H0 * self.E(z)
    
    def Comoving_Distance(self, z):
        """
        calculate the comoving distance
        z: redshifts
        """
        return self.cd_interp(z)

    def dVdz(self,z):
        """
        calculate the diffrential comoving volume dV/dz in units of Mpc^3
        z: redshifts
        """
        drdz = self.c / self.H(z)
        cd2 = self.Comoving_Distance(z)
        cd2 = cd2 * cd2
        dcv = 4 * np.pi * cd2 * drdz / (self.Mpc2cm * self.Mpc2cm * self.Mpc2cm)
        return dcv
    
    def Luminosity_Distance(self, z):
        """
        calculate the luminosity distance in Mpc
        z: redshifts
        """
        dl = (1 + z) * self.Comoving_Distance(z)
        return dl

    def Luminosity(self, z, f=1., dnu=1000.):
        """
        calculate the intrinsic luminosity from flux
        f: flux in unit of Jy
        dnu: intrinsic spectral width
        z: redshifts
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        lum = f * self.Jy2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return lum

    def Energy(self, z, flu=1.0, dnu=1000.):
        """
        calculate the intrinsic energy from fluence
        flu: fluence in units of Jy ms
        dnu: intrinsic spectral width, default at 1000 MHz
        z: redshift
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        ener = flu / (1+z) * self.Jyms2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return ener
    
    def LuminosityDistance_to_z(self, ld):
        """
        converting luminosity distance to redshift
        """
        return self.Ld2z(ld)
    
    def Luminosity_Distance_dimless(self, z):
        """
        calculate the dimension less luminosity distance
        z: redshifts
        """
        dl = (1 + z) * self.Comoving_Distance(z) / (2. * self.Cd1)
        return dl
    
    def DispersionMeasure_IGM(self, z, chi=7./8):
        """
        calculate dispersion measure of intergalactic medium by integrating redshift
        """
        return self.dmigm(z) * self.f_IGM * chi / self.Mpc2cm * 1e6
        
    def Luminosity_to_Flux(self, z, lum, dnu=1000):
        """
        calculate the flux density observed at a given luminosity
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld*ld
        flux = lum/4/np.pi/ld2/dnu/self.MHz2Hz/self.Jy2CGS
        return flux
    
    def Energy_to_Flu(self, z, ener, dnu=1000):
        """
        calculate the observed fluence when knowing intrinsic energy
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld*ld
        flu = ener*(1+z)/4/np.pi/ld2/dnu/self.MHz2Hz/self.Jyms2CGS
        return flu

    def DMeq(self, z, dme, dmhost):
        """
        The DM equation which is solved to get redshift value
        """
        dmi = self.DispersionMeasure_IGM(z)
        dmh = dmhost/(1+z)
        return dmi + dmh - dme

    def GetZ(self, dme, dmhost):
        """
        solving the differential equation to get redshift
        """
        z = fsolve(self.DMeq, self.z0, args=(dme, dmhost))
        z = float(z)
        return z

class AstroDistribution:
    def __init__(self):
        self.cos = Cosmology()
        self.vpar_etg = np.array([0.001713, 1.099, 0.2965, 0.01246, 1.055, 0.7262])
        self.vpar_ltg_ne2001 = np.array([0.01715, 1.062, 0.5202, 0.00416, 0.7227, 1.151])
        self.vpar_ltg_ymw16 = np.array([0.01561, 0.759, 0.3013, 0.01889, 1.042, 0.5791])
        self.vpar_alg_ne2001 = np.array([0.005485, 0.8665, 1.009, 0.01406, 1.069, 0.5069])
        self.vpar_alg_ymw16 = np.array([0.01199, 0.7597, 0.3082, 0.01735, 1.048, 0.6025])
        self.Zmax = 4
        self.Zmin = 2e-6
        self.DMsmax = 50.
        self.Epsmin = 0.5

    def Schechter_log(self, logl, phis, logls, alpha):
        """
        The Schecheter luminosity function per logarithmic luminosity
        logl: the logarithmic luminosity
        phis: normalization constant
        logls: the cut-off luminosity
        alpha: power index
        """
        l = np.power(10., logl)
        ls = np.power(10., logls)
        phi = np.log(10) * phis * np.power(l / ls, (alpha + 1)) * np.exp(-l / ls)
        return phi
    
    def log_Schechter_log(self, logl, logls, alpha, logl0):
        """
        The logarithm of Schecheter luminosity function per logarithmic luminosity
        logl: the logarithmic luminosity
        phis: normalization constant
        logls: upper cut-off luminosity
        alpha: power index
        logl0: lower cut-off luminosity
        """
        phi = (logl - logls) * (alpha+1) * np.log(10.) - np.power(10., logl-logls)
        lik = phi.copy()
        #lik [logl >= logl0] = (logl - logls) * (alpha+1) * np.log(10.) - np.power(10., logl-logls)
        lik [logl < logl0] = -1e99
        return lik

    def log_IntBeam(self, logl, logls, alpha, logl0):
        ratio = np.power(10., logl-logls)
        lik0 = gammainc(alpha+1, ratio) - gammainc(alpha+1, 2*ratio)
        lik = lik0/np.log(2)      # Beam efficiency from 50% to 100%
        ind = lik <= 0
        lik[ind] = 1e-199
        loglik = np.log(lik)
        #print loglik
        #loglik = np.ones(lik.shape) * -1e99
        #loglik[lik<=0] = -1e99
        loglik[logl < logl0] = -1e99
        #ind = lik>0 or logl>=logl0
        #loglik[ind] = np.log(lik[ind])
        return loglik

    def IntFluxBeam(self, logls, alpha, logl0):
        lr = np.power(10., logl0-logls)
        p1 = -1 + np.power(2., alpha+1)
        p2 = np.power(lr, alpha+1)
        p3 = spf.gamma(alpha+1)
        p4 = spf.gamma(alpha+2)
        lik = -p1*p2*p3*p3/p4/p4 + p3*np.log(2)
        return lik

    def IntLum(self, eps, logls, alpha, logl0):
        ratio = np.power(10., logl0-logls)
        return gammainc(alpha+1, ratio/eps)
    
    def Distribution_Local_galaxy_DM(self, dmv, vpar):
        """
        General form of DM distribution function of host galaxies, using double gaussian function in logarithmic DM
        vpar[i] is the i-th parameter of this distribution function
        """
        val = vpar[0] * np.exp(-np.power((dmv - vpar[1]) / vpar[2], 2.)) \
              + vpar[3] * np.exp(-np.power((dmv - vpar[4]) / vpar[5], 2.))
        return val

    def ThetaFunc(self, x):
        """
        Normalize the values positive
        """
        return 0.5 * (np.sign(x) + 1)
            
    def func_gaussian(self, dmv, vpar):
        """
        default gaussina distribution for non-galaxy case
        """
        dmoff = dmv - vpar[0]
        sig = vpar[1]
        sig = sig * sig
        return np.exp(-0.5 * dmoff * dmoff / sig) * self.ThetaFunc(dmv)
    
    def func_uniform(self, dmv, vpar):
        pdf = np.ones(dmv.shape)
        pdf[dmv >= vpar[1]] = 0
        pdf[dmv <= vpar[0]] = 0
        return pdf
    
    def Distribution_HostGalaxyDM(self, dmv0, fgalaxy_type=None, vpar=np.array([0,50])):
        """
        DM distribution functions of different type of host galaxies
        ETG: early-tyep galaxies
        LTG: late-type galaxies
        ALG: all the galaxies
        NE2001: the referenced galaxy electron density using NE2001 model
        YMW16: the referenced galaxyt electron density using YMW16 model
        """
        ind = dmv0 < 0
        if not fgalaxy_type: 
            fgalaxy_type = self.func_gaussian
        elif fgalaxy_type == 'ETG':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_etg)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'LTG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), 
                    self.vpar_ltg_ne2001)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'LTG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), 
                    self.vpar_ltg_ymw16)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'ALG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ne2001)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'ALG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ymw16)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        else:
            return fgalaxy_type(dmv0,vpar)
    
    
    def log_Distribution_HostGalaxyDM(self, dmv, fgalaxy_type=None, vpar=np.array([0,50])):
        if not fgalaxy_type:
            dmoff = dmv[dmv > 0] - vpar[0]
            sig = vpar[1]
            sig = sig * sig
            res = dmv.copy()
            res[dmv < 0] = -1e99
            res[dmv > 0] = -0.5 * dmoff * dmoff / sig
            return res
        else:
            val = self.Distribution_HostGalaxyDM(dmv, fgalaxy_type=fgalaxy_type)
            ind=val <= 0
            indv=val>0
            val[indv] = np.log(val[indv])
            val[ind] = np.ones(val[ind].shape) * -1e99
            return val
    
    def SFR(self, z):
        """
        Star-forming history, the values taken from Hopkins & Beacom (2016)
        """
        sfr = (0.017 + 0.13 * z)/(1 + np.power(z/3.3, 5.3))
        return sfr 

    def kappa(self, z):
        """
        Normalized SFH from redshift of z to redshift of 0 (nearby universe)
        """
        return np.sqrt(self.SFR(0)/self.SFR(z))

    def Distribution_volume(self, z):
        """
        Differential comoving volume
        """
        r = self.cos.Comoving_Distance(z)/self.cos.Comoving_Distance(1)
        pv = r * r / self.cos.E(z)
        return pv
    
    def log_Distribution_volume(self, z):
        """
        Logarithimc differential comoving volume
        """
        if type(z)==np.ndarray:
            ind=z<0
            pv=np.log(self.Distribution_volume(z))
            pv[ind]=-1e99
            return pv
        else:
            if z<0:
                return -1e99
            else:
                return np.log(self.Distribution_volume(z))

    def IntDMsrc(self, u1, u2, vpar):
        """
        Analytic integral when marginalizing distribution of DMsrc
        """
        a1 = vpar[0]
        b1 = vpar[1]
        c1 = vpar[2]
        a2 = vpar[3]
        b2 = vpar[4]
        c2 = vpar[5]
        k1 = np.power(10.,b1)*a1*c1*np.exp(c1*c1*np.log(10)*np.log(10)/4)
        k2 = np.power(10.,b2)*a2*c2*np.exp(c2*c2*np.log(10)*np.log(10)/4)
        q1 = (c1*c1*np.log(10)*np.log(10)+b1*np.log(100)-2*np.log(u1))/c1/np.log(100)
        q2 = (c1*c1*np.log(10)*np.log(10)+b1*np.log(100)-2*np.log(u2))/c1/np.log(100)
        q3 = (c2*c2*np.log(10)*np.log(10)+b2*np.log(100)-2*np.log(u1))/c2/np.log(100)
        q4 = (c2*c2*np.log(10)*np.log(10)+b2*np.log(100)-2*np.log(u2))/c2/np.log(100)
        int_h = np.log(10)*np.sqrt(np.pi)/2 * (k1*(spf.erf(q2)-spf.erf(q1))+k2*(spf.erf(q4)-spf.erf(q3)))
        int_hs = int_h/self.DMsmax   #   integral including uniform DMsrc
        return int_hs

    def log_IntDMsrc(self, u1, u2, gtype=None):
        """
        Logarithmic integrals of above maginalization in different galaxy type cases
        """
        if not gtype:
            gtype = self.func_gaussian
        elif gtype == 'ETG':
            res = np.log(self.IntDMsrc(u1, u2, self.vpar_etg))
            return res
        elif gtype == 'LTG_NE2001':
            res = np.log(self.IntDMsrc(u1, u2, self.vpar_ltg_ne2001))
            return res
        elif gtype == 'LTG_YMW16':
            res = np.log(self.IntDMsrc(u1, u2, self.vpar_ltg_ymw16))
            return res
        elif gtype == 'ALG_NE2001':
            res = np.log(self.IntDMsrc(u1, u2, self.vpar_alg_ne2001))
            return res
        elif gtype == 'ALG_YMW16':
            res = np.log(self.IntDMsrc(u1, u2, self.vpar_alg_ymw16))
            return res
        else:
            return -1e99
           
    def log_distr_fdmz(self, dnu, logflux, dme, z, logls, alpha, logl0, gtype=None):
        """
        Logarithmic joint distribution function of flux, DM and redshift
        """
        flux = np.power(10., logflux)
        logl = np.log10(self.cos.Luminosity(z, f=flux, dnu=dnu))
        logint1 = self.log_IntBeam(logl, logls, alpha, logl0)
        logfz = self.log_Distribution_volume(z)
        dmi = self.cos.DispersionMeasure_IGM(z)
        u1 = (dme-dmi)*(1+z)*self.kappa(z)
        u2 = ((dme-dmi)*(1+z)-self.DMsmax)*self.kappa(z)
        logint2 = np.ones(u2.shape) * (-1e99)
        ind = u2 > 0
        logint2[ind] = self.log_IntDMsrc(u1[ind],u2[ind],gtype=gtype)
        loglikv = logint1 + logfz + logint2 + np.log(1+z)
        return loglikv

    def log_distr_fdm(self, dnu, logflux, dme, logls, alpha, logl0, gtype=None):
        """
        Logarithmic joint distribution function of flux and DM after maginalization of redshift
        :param dnu: intrinsic bandwidth, i.e. 1GHz
        :param logflux: logarithmic flux
        :param dme: dme
        :param logls: logarithmic cut-off luminosity
        :param alpha: lf index
        :param logl0: minimum of lf L
        :param fgalaxy_type: galaxy type
        :return logarithmic likelihood:
        """
        #stepdms = 100/1000.
        #vdms = np.arange(0, 100, stepdm)
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        lik = 0
        for z in vz:
            likv = np.exp(self.log_distr_fdmz(dnu, logflux, dme, z, logls, alpha, logl0, gtype=gtype))
            lik += z * stepz * likv
        ind = lik > 0
        ind2 = lik <= 0
        loglik = lik.copy()
        loglik[ind] = np.log(lik[ind])
        loglik[ind2] = np.ones(loglik[ind2].shape) * -1e99
        return loglik
    
    def Norm1D(self, logflux_thre, dnu, logls, alpha, logl0):
        """
        Normalization factor
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        stepeps = (1-0.5)/500.
        veps = np.arange(0.5, 1, stepeps)
        nf = 0
        #Lthre = self.cos.Luminosity(vz, f=np.power(10.,logflux_thre), dnu=dnu)
        #loglt = np.log10(Lthre)
        #ind = loglt < logl0
        #loglt[ind] = logl0
        for z in vz:
            Lthre = self.cos.Luminosity(z, f=np.power(10.,logflux_thre), dnu=dnu)
            loglt = np.log10(Lthre)
            if loglt < logl0:
               loglt = logl0
            lika = self.Distribution_volume(z)
            likb = self.IntLum(veps, logls, alpha, loglt*np.ones(veps.shape))/veps/np.log(2)
            nf += z*stepz*lika*np.sum(stepeps*likb)
        #part_a = np.exp(self.log_Distribution_volume(vz))
        #part_a = self.Distribution_volume(vz)
        #part_b = gammainc(alpha+1, np.power(10., loglt-logls))
        #print part_b
        #lik = part_a*part_b
        #print lik
        #nf = np.sum(vz*lik*stepz)
        #print nf
        #lognf = np.log(nf)
        return nf

class Loadfiles:
    def __init__(self):
        self.MHz2Hz = 1e6
        self.ms2s = 1e-3

    def LoadCatalogue(self, fname):
        cat = np.loadtxt(fname, dtype='string')
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]
        
        cat2['Tel'] = np.array(cat2['Tel'])
        #cat2['STD'] = np.array(cat2['STD'], dtype=float)
        cat2['SEFD'] = np.array(cat2['SEFD'], dtype=float)
        cat2['Npol'] = np.array(cat2['Npol'], dtype=float)
        cat2['SIG'] = np.array(cat2['SIG'], dtype=float)
        cat2['BW'] = np.array(cat2['BW'], dtype=float)
        cat2['DME_NE2001'] = np.array(cat2['DME_NE2001'], dtype=float)
        cat2['DME_YMW16'] = np.array(cat2['DME_YMW16'], dtype=float)
        cat2['W'] = np.array(cat2['W'], dtype=float)
        cat2['S'] = np.array(cat2['S'], dtype=float)
        cat2['Seu'] = np.array(cat2['Seu'], dtype=float)
        cat2['Sel'] = np.array(cat2['Sel'], dtype=float)
        cat2['F'] = np.array(cat2['F'], dtype=float)
        cat2['Feu'] = np.array(cat2['Feu'], dtype=float)
        cat2['Fel'] = np.array(cat2['Fel'], dtype=float)
        return cat2

    def SEFD(self, g, tsys):
        sefd = tsys/g
        return sefd

    def RMEq(self, snr, sefd, npol, bw, w):
        s = snr*sefd/np.sqrt(npol*bw*self.MHz2Hz*w*self.ms2s)
        return s

    def LoadSimuData(self, fname):
        cat = np.loadtxt(fname, dtype='string',comments='XXX')
        row, col = cat.shape
        cat2 = {}
        cat[0,0]=cat[0,0][1:]
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]
        
        cat2['DMe'] = np.array(cat2['DMe'], dtype=float)
        cat2['S'] = np.array(cat2['S'], dtype=float)
        cat2['DMi'] = np.array(cat2['DMi'], dtype=float)
        cat2['DMh'] = np.array(cat2['DMh'], dtype=float)
        cat2['DMs'] = np.array(cat2['DMs'], dtype=float)
        cat2['Z'] = np.array(cat2['Z'], dtype=float)
        cat2['logL'] = np.array(cat2['logL'], dtype=float)
        return cat2

def gammainc(alpha, x):
    if alpha==0:
        return -spf.expi(-x)

    elif (alpha<0):
        return (gammainc(alpha+1,x)-np.power(x, alpha)*np.exp(-x))/alpha

    else:
        return spf.gammaincc(alpha,x)*spf.gamma(alpha)

def getargv(argv, key):
    for i in range(0, len(argv)):
        arg = argv[i]
        if (arg == key):
            return argv[i + 1]

def chkargv(argv, key):
    for i in range(0, len(argv)):
        arg = argv[i]
        if (arg == key):
            return True
    return False

def Sampling1D(x, y, x1, x2, n):
    nt = 0
    res = np.array([])
    while (nt<n):
        fuc = interp1d(x, y / np.max(y))
        vx = np.random.uniform(x1, x2, n-nt)
        vy = np.random.uniform(0, 1, n-nt)
        res = np.append(res, vx[vy <= fuc(vx)])
        nt = len(res)
    return res

def SamplingND(fuc, par_range, maxv_ori, n):
    nt = 0
    maxv = maxv_ori
    res = np.array([])
    npar, m = par_range.shape
    res = res.reshape((0, npar))
    while (nt < n):
        vpar = np.random.uniform(0, 1, (n-nt, npar))
        for i in range(npar):
            lv = par_range[i,0]
            rv = par_range[i,1]
            vpar[:,i] = vpar[:,i] * (rv-lv) + lv
        
        vy = np.random.uniform(0, 1, n-nt)*maxv
        fv = fuc(vpar)
        if (np.max(fv) > maxv):
            maxv = np.max(fv)*1.5
            nt = 0
            res = np.array([])
            res = res.reshape((0, npar))
        else:
            res = np.vstack( (res, vpar[vy <= fv,:]))
            nt, m = res.shape
        print nt
    return res


def itp2dfunc(x, y, mx, my, mz):
    kx = mx[0,1] - mx[0,0]
    kx1 = (mx[0,1] - x)/kx
    kx2 = (x - mx[0,0])/kx
    ky = my[1,0] - my[0,0]
    ky1 = (my[1,0] - y)/ky
    ky2 = (y - my[0,0])/ky
    fxy1 = kx1*mz[0,0] + kx2*mz[1,0]
    fxy2 = kx1*mz[0,1] + kx2*mz[1,1]
    fxy = ky1*fxy1 + ky2*fxy2
    return fxy

def findgrid(x, y, mx, my, mz):
    for j in range(my.shape[0]):
        if y>=my[j,0] and y<=my[j+1,0]:
            ny = j
    for i in range(mx.shape[1]):
        if x>=mx[0,i] and x<=mx[0,i+1]:
            nx = i
    return [nx, ny]

def itp2d(x, y, mx, my, mz):
    mxnew, mynew = np.meshgrid(x,y)
    res = np.ones(mxnew.shape)
    if np.min(mxnew)>=np.min(mx) and np.max(mxnew)<=np.max(mx) and np.min(mynew)>=np.min(my) and np.max(mynew)<=np.max(my):
        for i in range(mxnew.shape[0]):
            for j in range(mxnew.shape[1]):
                px = mxnew[i,j]
                py = mynew[i,j]
                nx = findgrid(px, py, mx, my, mz)[0]
                ny = findgrid(px, py, mx, my, mz)[1]
                mxx, myy=np.meshgrid([mx[0, nx], mx[0, nx+1]], [my[ny, 0], my[ny+1, 0]])
                mzz = mz[ny:ny+2, nx:nx+2]
                res[i,j] = itp2dfunc(px, py, mxx, myy, mzz)
        return res
    else:
        print "Error: The inputs are out of interpolation range."
