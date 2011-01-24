#!/usr/bin/env python
#
# Basic numerical implementation of the Mansoori equation for a binary HS fuild.
#
# REF: The Journal of Chemical Physics --  Volume 54, Issue 4, pp. 1523-1525 
#

import sys
from math import log,pi,fabs
from AlderCrystal import *

class Mansoori:
    """
    Calculates pressures and free energies (wrt IG) as a function of density for the
    Mansoori semi-empirical equation for a binary hard sphere fluid.
    ---
    Length unit is \sigma_A.
    """
    
    def __init__(self, alpha, xB):
        "Initialise with the desired radius ratio and fraction composition of B (ie xA = 1.0 - xB)"
        sA = 1.0
        sB = alpha*sA
        xA = 1.0 - xB
        self.alpha = alpha
        self.xA = xA
        self.xB = xB
        self.sA = sA
        self.sB = sB
        
        # Determine the parameters y1,y2,y3 from the above:
        self.y1 = (xA*xB*(sA+sB)*((sA-sB)**2))/(xA*(sA**3) + xB*(sB**3))
        self.y2 = (xA*xB*sA*sB*(xA*(sA**2)+xB*(sB**2))*((sA-sB)**2))/((xA*(sA**3) + xB*(sB**3))**2)
        self.y3 = ((xA*(sA**2)+xB*(sB**2))**3)/((xA*(sA**3) + xB*(sB**3))**2)
        
    def Z(self, vf):
        "Returns the compressibility Z=PV/NkT at the given volume fraction."
        # First compute the compressibility in the form \beta.P/\rho
        Z = (1.0 + vf + (vf**2.0) - 3.0*vf*(self.y1+vf*self.y2) - self.y3*(vf**3.0))/((1.0-vf)**3.0)
        return Z
    
    def V(self, vf):
        "Calculate the volume [1/N.sigma_A] from the volume fraction."
        # Two lines, like this one, are suspect!  Other one has been replaced ATM.
        return (pi/(6.0*vf))*(self.xA*(self.sA**3.)+self.xB*(self.sB**3.))
        #return (pi/(6.0*vf))
    
    def p(self, vf):
        "Returns the pressure at the given volume fraction."
        # First compute the pressure in the form \beta.P/\rho
        bPod = self.Z(vf)
        # Now change to our units such that p = \beta.P.\sigma_A^3
        bPsA = bPod / self.V(vf)
        return bPsA
    
    def FEld(self,vf):
        "Return the Helmholtz free energy (wrt IG) at the given volume fraction."
        # Form from Eldridge AB13 paper - this does not behave well, and I can't derive it from Mansoori.
        fh = ( 3*(1-self.y1)/2 + 3*(0.5-vf)*self.y2 + (3.0/2.0 - 2.0*vf)*self.y3 )/((1-vf)**2) - log(1.-vf) + self.y3*log(1.-vf)
        return fh
    
    def FH(self,vf):
        "Return the Helmholtz free energy (wrt IG) at the given volume fraction."
        # Form from Eldridge AB13 paper - this does not behave well, and I can't derive it from Mansoori.
        #fh = ( 3.0*(1.0-self.y1)/2.0 + 3.0*(0.5-vf)*self.y2 + (3.0/2.0 - 2.0*vf)*self.y3 )/((1.0-vf)**2.) - log(1.0-vf) + self.y3*log(1.0-vf)
        # Form from original Mansoori paper.
        # This agrees with Carnahan-Starling equation at alpha=1.0, w_B = 0.0.
        fh = -3.*(1.-self.y1+self.y2+self.y3)/2. + (3.*self.y2+2.*self.y3)/(1.-vf) +  3.*(1.-self.y1-self.y2-self.y3/3.)/(2.*((1-vf)**2.)) + (self.y3-1.)*log(1.-vf)
        # This expression derived from Maple integration of E.O.S. Same as Mansoori form.
        #fh = 3.0*(1.0-self.y1-self.y2-self.y3/3.0)/((2.0*(vf-1)**2))
        #fh += (2.0*self.y3+3.0*self.y2)/(1.0-vf)
        #fh += (self.y3-1.0)*log(1.0-vf)
        #fh += 3.0*(self.y1-self.y2-self.y3-1.0)/2.0
        return fh
    
    def S(self,vf):
        "Return the excess entropy, in units of [Nk]."
        # The reason for this ln(Z) term is not clear - perhaps a choice og IG entropy.
        # This reproduces the results in the Mansoori tables.
        # This should not be used for my work.
        return -self.FH(vf)+log(self.Z(vf))
    
    def f(self,vf):
        "Return the excess Helmholtz free energy, in units of [NkT]."
        # Returning this as is gives agreement with the CS form.
        fe = self.FH(vf)
        if( self.xA > 0.0 ):
            fe = fe + self.xA*log(self.xA)
        if( self.xB > 0.0 ):
            fe = fe + self.xB*log(self.xB)
        return fe
    
    def g(self,vf):
        "Return the Gibbs free energy (w.r.t IG), in units of [NkT]."
        return self.f(vf) + self.Z(vf) - 1.0 - log(self.Z(vf)) 
    #- log(self.Z(vf)) -log(1.0/self.p(vf))
        #- log(self.Z(vf))

    def g_abs(self,vf):
        "Return the Gibbs free energy (absolute), in units of [NkT]."
        return self.f(vf) + self.Z(vf) - 1.0 + log(self.V(vf)) 
    
    def GH(self,vf):
        "Return the Gibbs free energy (wrt IG) at the given volume fraction. G = F + PV."
        # This agrees with the tabulated results in the Mansoori paper, but the ln(Z) term is in there.
        return self.Z(vf) - 1 - self.S(vf)
    
    def vf(self, pres):
        "Inverse operation on the EOS, find the density for a given pressure."
        # Use Newton-Raphson approach:
        vf = 0.7
        tol = 1.0e-6
        delta = 1.0e-7
        itry = 0
        max_try = 20
        while fabs(self.p(vf)-pres) > tol and itry < max_try:
            #print "Crystal; ",itry," Want ",pres," Got ",self.p(vf), " at vf=",vf
            vf = vf - (self.p(vf)-pres)/((self.p(vf+delta)-self.p(vf))/delta)
            itry += 1
            
        if itry >= max_try:
            print "WARNING: Mansoori vf root finder exited early after ",itry," iterations."
        # Return the vf:
        if vf < 0.0:
             return 0.0
        #if vf > 0.74:
        #    print "Got vf=",vf," resetting..."
        #    return 0.74
        return vf

    def uExpr(self,vf,t1,t2,t3,nuOx):
        "Returns the chemical potential based on the formula in Bartlett [1990]."
        u = 3.0*(t1-t2-t3-1)/2.0 + (t3-1.0)*log(1.0-vf) + (2.0*t3+3.0*t2)/(1.0-vf)
        u += (3.0*(1-t1-t2)-t3)/(2.0*((1.0-vf)**2.0))
        upart = (1.0-self.y3)/(1.0-vf) + (2.0*self.y3+3.0*self.y2)/((1.0-vf)**2.0)
        upart += (3.0*(1.0-self.y1-self.y2)-self.y3)/((1.0-vf)**3.0)
        u += nuOx*upart
        return u
    
    def uA(self,vf):
        "Return the chemical potential of the species A at the given volume fraction & B-conc."
        # Based on formulae from Bartlett [1990]
        # Get x1, x2 from xB
        x2 = self.xB
        x1 = 1.0-x2
        # And the sigmas
        s1 = self.sA
        s2 = self.sB
        # And the ys (already defined upon instanciation)
        #self.y1 = (x1*x2*(s1+s2)*((s1-s2)**2))/(x1*(s1**3) + x2*(s2**3))
        #self.y2 = (x1*x2*s1*s2*(x1*(s1**2)+x2*(s2**2))*((s1-s2)**2))/((x1*(s1**3) + x2*(s2**3))**2)
        #self.y3 = ((x1*(s1**2)+x2*(s2**2))**3)/((x1*(s1**3) + x2*(s2**3))**2)
        # And the nus (partial volume fractions)
        if x1 == 0:
            nuA = 0.0
            nuB = vf
            nuOx = 0.0
        elif x2 == 0:
            nuA = vf
            nuB = 0.0
            nuOx = nuA/x1
        else:
            #nuA = vf*(s1**3.)*x1/((s1**3)*x1+(s2**3.)*x2)
            #nuB = vf*(s2**3.)*x2/((s1**3)*x1+(s2**3.)*x2)
            nuA = vf/( 1.0 + (x2*(s2**3.)/(x1*(s1**3.))) )
            nuB = vf/( 1.0 + (x1*(s1**3.)/(x2*(s2**3.))) )
            nuOx = nuA/x1
        # Define the 'theta' parameters dy/dNA|NB
        sigs = ((s1-s2)**2.0)*(s1+s2)
        # By hand:
        #t1 = self.y1 + (sigs/(x1*(s1**3.)+x2*(s2**3.)))*(x2 - (x1*x2/(x1*(s1**3.)+x2*(s2**3.)))*(x1*(s1**3.)+x2*(s2**3.)+(s1**3.)))
        #t2 = self.y2 + (s1*s2*((s1-s2)**2)*x2/((x1*s1**3+x2*s2**3)**2))*(2*x1*(s1**2)+x2*(s2**2)-x1*(x1*(s1**2)+x2*(s2**2))-2*x1*(s1**3)*(x1*(s1**2)+x2*(s2**2))/(x1*(s1**3)+x2*(s2**3)))
        #t3 = self.y3 + ( ((x1*(s1**2)+x2*(s2**2))**2)/((x1*(s1**3)+x2*(s2**3))**2) )*(3*(s1**2)-(x1*(s1**2)+x2*(s2**2))*(1+2*(s1**3)/(x1*(s1**3)+x2*(s2**3)) ))
        # Via maple:
        t1 = self.y1 - sigs*x2*(x1**2*s1**3-x2**2*s2**3)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**2)
        t2 = self.y2 - (s1*s2*((s1-s2)**2))*x2*(x1**3*s1**5-x1**2*s1**2*x2*s2**3-2.*x1*s1**2*x2**2*s2**3+2.*x2*s2**2*x1**2*s1**3+x2**2*s2**2*x1*s1**3-x2**3.*s2**5)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**3)
        t3 = self.y3 + (x1*s1**2+x2*s2**2)**2*x2*(2*s1**2*x1*s2**3+s1**5*x1+3.*s1**2*x2*s2**3-3.*s2**2*x1*s1**3-x2*s2**5-2.*s1**3*x2*s2**2)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**3)
        # Determine this chemical potential:
        return log(self.xA/self.V(vf)) + self.uExpr(vf, t1, t2, t3, nuOx)
    
    def uB(self,vf):
        "Return the chemical potential of the species B at the given volume fraction & B-conc."
       # Based on formulae from Bartlett [1990]
        # Get x1, x2 from xB
        x2 = self.xB
        x1 = 1.0-x2
        # And the sigmas
        s1 = self.sA
        s2 = self.sB
        # And the ys (already defined upon instanciation)
        #self.y1 = (x1*x2*(s1+s2)*((s1-s2)**2))/(x1*(s1**3) + x2*(s2**3))
        #self.y2 = (x1*x2*s1*s2*(x1*(s1**2)+x2*(s2**2))*((s1-s2)**2))/((x1*(s1**3) + x2*(s2**3))**2)
        #self.y3 = ((x1*(s1**2)+x2*(s2**2))**3)/((x1*(s1**3) + x2*(s2**3))**2)
        # And the nus (partial volume fractions)
        if x1 == 0.0:
            nuA = 0.0
            nuB = vf
            nuOx = nuB/x2
        elif x2 == 0.0:
            nuA = vf
            nuB = 0.0
            nuOx = 0.0
        else:
            #nuA = vf*(s1**3.)*x1/((s1**3.)*x1+(s2**3.)*x2)
            #nuB = vf*(s2**3.)*x2/((s1**3.)*x1+(s2**3.)*x2)
            nuA = vf/( 1.0 + (x2*(s2**3)/(x1*(s1**3))) )
            nuB = vf/( 1.0 + (x1*(s1**3)/(x2*(s2**3))) )
            nuOx = nuB/x2
        # Define the 'theta' parameters dy/dNA|NB
        sigs = ((s1-s2)**2.0)*(s1+s2)
        # By hand:
        # Via maple:
        t1 = self.y1 + sigs*x1*(x1**2*s1**3-x2**2*s2**3)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**2)
        t2 = self.y2 + (s1*s2*((s1-s2)**2))*x1*(x1**3*s1**5-x1**2*s1**2*x2*s2**3-2*x1*s1**2*x2**2*s2**3+2*x2*s2**2*x1**2*s1**3+x2**2*s2**2*x1*s1**3-x2**3*s2**5)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**3)
        t3 = self.y3 - (x1*s1**2+x2*s2**2)**2*x1*(-3*s2**2*x1*s1**3-s2**5*x2-2*s2**2*x2*s1**3+x1*s1**5+3*s1**2*x2*s2**3+2*s2**3*x1*s1**2)/((x1+x2)**2*(x1*s1**3+x2*s2**3)**3)
        # Determine this chemical potential:
        return log(self.xB/self.V(vf)) + self.uExpr(vf, t1, t2, t3, nuOx)
    
    def xB_for_coex_pres(self, pres, Gcry, numB):
        # For the given number of B particles for each A (numB) and the overall Gibbs free energy per particle Gcry:
        # Determine the value of xB that would lead to coexistance at the given pressure.
        # Use coarse linear cross-over-spotter to find rough value
        Xr = 1.0e-10
        step = 0
        steps = 10000
        u_tox = 0.0
        Xr_u = 0.0
        while step < steps:
            man = Mansoori(self.alpha,Xr)
            man_vf = man.vf(pres)
            man_utot = man.uA(man_vf) + numB*man.uB(man_vf)
            u_diff = man_utot - (1.0+numB)*Gcry
#            print "UCALC: ",pres," ",Xr," ",man.uA(man_vf)," ",man.uB(man_vf)," ",Gcry," ",man_utot," ",(1.0+numB)*Gcry, " ", u_diff
            if step == 0:
                u_tox = u_diff
                Xr_u = Xr
            else:
                if u_tox/u_diff < 0.0:
                    Xr_u = Xr
                    u_tox = u_diff
            # Next:
            step = step + 1
            Xr = (1.0-(1.0e-10))*step/steps

#        print "(init: ", pres, " ", u_tox, " ", u_diff, " guess: ", Xr_u
        #return Xr_u
        
        # Use Newton-Raphson approach to refine:
        Xr = Xr_u
        tol = 1.0e-6
        delta = 1.0e-7
        itry = 0
        max_try = 10000
        u_diff = 2.0*tol
        while fabs(u_diff) > tol and itry < max_try:
            # Construct Mansoori fluid at the given composition
            man = Mansoori(self.alpha,Xr)
            man_vf = man.vf(pres)
            man_utot = man.uA(man_vf) + numB*man.uB(man_vf)
            u_diff = man_utot - (1.0+numB)*Gcry
            # Now set up results at +delta:
            manDelta = Mansoori(self.alpha,Xr+delta)
            manDelta_vf = manDelta.vf(pres)
            manDelta_utot = manDelta.uA(manDelta_vf) + numB*manDelta.uB(manDelta_vf)
            uDelta_diff = manDelta_utot - (1.0+numB)*Gcry
            # Now apply N-R:
            #print pres, " ", Xr, " ", man_utot," ", manDelta_utot, " ", u_diff, " ", (u_diff)/((manDelta_utot-(1.0+numB)*Gcry)/delta)
            Xr = Xr - delta*u_diff/uDelta_diff
            if Xr <= 0.0:
                Xr = 0.0+tol/100.0
            if Xr >= 1.0:
                Xr = 1.0-tol/100.0
                
            itry += 1
            
        if itry >= max_try:
            print "# WARNING: Mansoori xB root finder exited early after ",itry," iterations."
            
        # Return the vf:
        return Xr
    
# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) >= 2:
        alpha = float(sys.argv[1])
    else:
        alpha = 1.0
    # Choose the fraction that is B:
    if len(sys.argv) >= 3:
        nB = float(sys.argv[2])
        if nB < 0.0:
            xB = - nB
        else:
            xB = nB/(1+nB)
            
        max_xB = xB
        
    else:
        # Scan the whole range of xB
        xB = 0.0
        max_xB = 1.0
        
    steps = 501
    max_vf = 0.65
    print "# {@alpha=",alpha,"} vf x_B PV/NkT P [s_A/kT] f [Nkt] g[NkT] g_abs [NkT]"
    cryA = AlderCrystal(alpha)
    while xB < max_xB+0.0000001:
        man = Mansoori(alpha, xB)
        # Loop over the density:
        step = 1
        while step <= steps:
            vf = max_vf*step/steps
            pres = man.p(vf)
            vf_cry = cryA.vf(pres)
            Gcry = cryA.g_abs(vf_cry)
            if pres > 12:
#                print "(G ",pres, " " ,Gcry, " ", 
                Xr_t = man.xB_for_coex_pres(pres, Gcry, 0.0)
                xB_c = 1.0 - (xB/Xr_t)
                xB_f = 1.0 - xB_c
                
            if pres <= 10 or Xr_t < xB:
                Xr_t = xB
                xB_c = 0.0
                xB_f = 1.0
            
            man_t = Mansoori(alpha, Xr_t)            
            vf_t = man_t.vf(pres)
            
            print vf, " ", xB, " ", man.Z(vf), " ", pres, " ", man.f(vf), " ", man.g(vf), " ", man.g_abs(vf), " ", man_t.g_abs(vf), " ", (1.0-xB)*Gcry, " ", Xr_t, " ", xB_c, " ", (xB_c*Gcry + xB_f*man_t.g_abs(vf_t)), " ", (xB_c*cryA.g(vf_cry) + xB_f*man_t.g(vf_t))
            step += 1
            
        print " "
        xB += 1.0/steps
        