#!/usr/bin/env python
#
# Basic numerical implementation of the Alder equation for a pure HS crystal.
#
# REF: The Journal of Chemical Physics -- Volume 70, Issue 1, pp. 473-481
#

import sys
from math import log,sqrt,pi,e,fabs
from Mansoori import *
from CSFluid import *

class AlderCrystal:
    """
    Calculates pressures and free energies (wrt IG) as a function of density for the
    Alder et alsemi-empirical equation for a hard sphere crystal.
    
    Changes from 'configurational Helmholtz free energy to w.r.t IG by adding one:
          (F-F_ig)/NkT = F/NkT + 1
          This is correct because (ignoring the internal energy) the free energy IG is:
           F_ig = -NkT
    ---
    Length unit is \sigma_A.
    """
    
    def __init__(self, alpha=1.0, xB=0.0):
        "Initialise parameters for the single crystal."
        # Store the crystal unit-scaling parameters:
        if alpha == 1.0:
            # If alpha is 1.0, then we are a pure A crystal:
            self.xA = 1.0
            self.sA = 1.0
            self.xB = 0.0
            self.sB = 1.0
            self.alpha = 1.0
        else:
            # Otherwise, this is a pure B crystal with reduced radius alpha:
            self.xA = 1.0-xB
            self.sA = 1.0
            self.xB = xB
            self.sB = alpha
            self.alpha = alpha
            
        
    def Z(self, vf):
        "Returns PV/NkT at the given volume fraction."
        # First compute the pressure in the form \beta.P/\rho
        v = self.Vred(vf)
        alv = v-1.0
        bPod = (3.0/alv) + 2.566 + 0.55*alv - 1.19*(alv*alv) + 5.95*(alv*alv*alv)
        return bPod
    
    def V(self, vf):
        "Calculate the volume [1/N.sigma_A] for this volume fraction."
        Vr = pi/(6.0*vf)
        return Vr
        
    
    def Vred(self,vf):
        "Calculate the reduced volume used in the original paper."
        return sqrt(2.0)*self.V(vf)
    
    def p(self, vf):
        "Returns the pressure at the given volume fraction."
        # First compute the pressure in the form \beta.P/\rho
        bPod = self.Z(vf)
        # Now change to our units, via \rho.\sigma_A^3, such that p = \beta.P.\sigma_A^3
        Vol = self.V(vf)
        Vol = self.V(vf)*((1.0-self.xB)+self.xB*(self.alpha**3))
        bPsA = bPod / Vol
        return bPsA
    
    def f(self,vf):
        "Return the Helmholtz free energy (wrt IG) at the given volume fraction."
        return self.fM(vf)
        #v = self.Vred(vf)
        #if v <= 1.0:
        #    v = 1.00000000001
        #return -3.0*log((v-1.0)/v) + 5.124*log(v) - 20.78*v + 19.04*(v**2.0)/2.0 -5.95*(v**3.0)/3.0 + 15.05 + 1
    
    def fM(self,vf):
        "Return the Helmholtz free energy (wrt IG) at the given volume fraction."
        v = self.Vred(vf)
        if v <= 1.0:
            v = 1.00000000001
        # This is the form in Young & Alder JChemPhys (70) p473:
        #return -3.0*log((v-1.0)/v) + 5.124*log(v) - 20.78*v + 9.52*(v**2.0) -1.98*(v**3.0) + 15.05 + 1
        # Add minor correction to constant to get sensible transition point against the Mansoori/CS fluid.
        #return -3.0*log((v-1.0)/v) + 5.124*log(v) - 20.78*v + 19.04*(v**2.0)/2.0 -5.95*(v**3.0)/3.0 + 15.05 + 1
        return -3.0*log((v-1.0)/v) + 5.124*log(v) - 20.78*v + 19.04*(v**2.0)/2.0 -5.95*(v**3.0)/3.0 + 15.05 - 0.075 + 1
        # This form gives sane results like those for AB13 by Eldridge/Frenkel etc.
        #return -3.0*log((v-1.0)/v) + 5.124*log(v) - 20.78*v + 19.04*(v**2.0)/2.0 -5.95*(v**3.0)/3.0 + 15.05 + 0.0564 + 1
        # From Maple - same result?:
        #return 1.9833333333333333333333333*(v**3)-9.52*(v**2)+20.78*v+3.*log(v-1.)-8.124*log(v) - 1
    
    def g(self,vf):
        "Return the Gibbs free energy (wrt IG) at the given volume fraction. G = F + PV."
        return self.gM(vf)
        #PV = self.Z(vf)
        #helmf = self.fM(vf)
        #return helmf + PV - 1

    def gM(self,vf):
        "Return the Gibbs free energy (wrt IG) at the given volume fraction. G = F + PV."
        # Include the correction to fit with the Mansoori/CS form:
        return self.Z(vf) + self.fM(vf) - 1.0 - log(self.Z(vf))
    #- log(self.Z(vf)) - log(1.0/self.p(vf))
        #- log(self.Z(vf)), -0.2

    def g_abs(self,vf):
        "Return the Gibbs free energy (absolute) at the given volume fraction. G = F + PV."
        # Include the correction to fit with the Mansoori/CS form:
        return self.Z(vf) + self.fM(vf) - 1.0 + log(self.V(vf))

    def f_fluid(self, vf):
        "Return the free-energy (wrt IG) of the Alder HS fluid at the given volume fraction."
        # This does not appear to be correct.  Give very wrong vf = 0 behaviour (-infinity).
        return (4*vf-3*(vf**2))/((1-vf)**2)+log(vf)+log(6.0/(pi*e)) + 1

    def p_fluid(self,vf):
        "Return the pressure of the Carnahan-Starling Fluid."
        csf = CSFluid(self.alpha)
        return csf.p(vf)
    
    def g_fluid(self,vf):
        "Return the Gibbs free energy (wrt IG) of the Alder fluid at the given volume fraction. G = F + PV."
        csf = CSFluid(self.alpha)
        PV = csf.Z(vf)
        helmf = self.f_fluid(vf)
        return helmf + PV - 1
    
    def vf(self, pres):
        "Inverse operation on the EOS, find the density for a given pressure."
        # Use simple iterative approach:
        vf = 0.45
        vfo = vf
        # If the requested pressure is too low, exit:
        if self.p(vf) > pres:
            return vf
        # Go up in pressure and vf, approaching the limit carefully:
        tries = 0
        res = 1000
        while self.p(vf) < pres and tries < 10000:
            vfo = vf
            #print vfo, " ", pres, " ", self.p(vfo)
            vf += (0.740-vf)/res
            tries+=1
        # Simplest thing is to return this crossover point - not even an interpolation:
        #return vf
        # Otherwise, interpolate linearly from vfo to vf:
        return vfo + (vf-vfo)*(pres-self.p(vfo))/(self.p(vf)-self.p(vfo))
    
    def u(self,vf):
        "Chemical potential for the single-phase crystal."
        # Expression taken from Bartlett [1990].
        v = self.Vred(vf)
        u = -3*log(v-1) + 9.124*log(v) - 9.52*v*v + 3.966*v*v*v +3*v/(v-1) + 6.5794 
        # This may need tweaking by + 0.103 to get chempots to agree in Tester.py
        return u
                
# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) >= 2:
        alpha = float(sys.argv[1])
    else:
        alpha = 1.0
    if len(sys.argv) == 3:
        nB = float(sys.argv[2])
        xB = nB/(1+nB)
    else:
        xB = 1.0
    xA = 1.0 - xB
    print "#Got alpha=",alpha," x_B=",xB, " x_A=",xA
    
    manA = Mansoori(1.0,0.0)
    manB = Mansoori(alpha,1.0)
    cryA = AlderCrystal(1.0, 0.0)
    cryB = AlderCrystal(alpha, 1.0)
    cry = AlderCrystal(alpha,xB)
    csf = CSFluid(alpha)
    pres_cx = 8.15
    # Loop over the density:
    steps = 500
    max_vf = 0.7405
    step = 1
    print "# {alpha=",alpha,"}"
    #," - vf, p, f, g, f_fluid, g_fluid"
    while step < steps:
        vf = 0.42+(max_vf-0.42)*step/steps
        pres = manA.p(vf)
        if pres >= pres_cx*sqrt(2.0):
            vfa = cryA.vf(pres)
            Ag = cryA.g(vfa)
            Af = cryA.f(vfa)
        else:
            vfa = manA.vf(pres)
            Ag = manA.g(vfa)
            Af = manA.f(vfa)
              
        if pres >= pres_cx*sqrt(2.0)/(alpha*alpha*alpha):
            vfb = cryB.vf(pres)
            Bg = cryB.g(vfb)
            Bf = cryB.f(vfb)
        else:
            vfb = manB.vf(pres)
            Bg = manB.g(vfb)
            Bf = manB.f(vfb)

        print vf, " ", manA.Z(vf), " ", pres, " ", (xA*Af + xB*Bf), " ", (xA*Ag + xB*Bg), " ", (cryA.f_fluid(vf)), " ", cryB.g_fluid(vf), " ", cry.Z(vf), " ", cry.p(vf), " ", (cry.Z(vf)/cry.Vred(vf)), " ", cry.g(vf), " ", csf.p(vf), " ", csf.g(vf), " ", vfa, " ", cryA.f(vfa)
#        print vf, cryA.V(vf)," ",cryA.Vred(vf)," ",cryA.Z(vf)," ",cryA.p(vf)
        step += 1
