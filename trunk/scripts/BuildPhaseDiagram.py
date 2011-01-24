#!/usr/bin/env python
#
# Pulls the literature equations and LS data together to produce phase diagrams and other overall results.
#
#

import sys
import math as m
from CSFluid import *
from AlderCrystal import *
from Mansoori import *

# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) >= 2:
        alpha = float(sys.argv[1])
    else:
        print "I need at least one argument! <Alpha> [<#Bs for each A>]"
        sys.exit()
        
    if len(sys.argv) == 3:
        nB = float(sys.argv[2])
        xB = nB/(1.0+nB)
    else:
        xB = 1.0
    xA = 1.0 - xB
    
    print "#Got alpha=",alpha," x_B=",xB, " x_A=",xA

    # End of CLARG parsing.

    # Produce GNU-plottable 3-D data for the different phases at the chosen alpha.
    # Data here is x_B, pres, Gibbs
    steps = 200
    
    # Mansoori/CS Fluid A    
    manA = Mansoori(1.0, 0.0)
    # Loop over the density:
    max_vf = 0.7405
    step = 1
    # Store EOS data:
    manA_vf = {}
    manA_p = {}
    while step < steps:
        vf = max_vf*step/steps
        #print >>f, 0.0, " ", cryA.p(vf), " ", cryA.g(vf)
        manA_vf[step] = vf
        manA_p[step] = manA.p(vf)
        step += 1
    # Mansoori/CS Fluid B
    manB = Mansoori(alpha, 1.0)
    # Loop over the density:
    max_vf = 0.7405
    step = 1
    # Store EOS data:
    manB_vf = {}
    manB_p = {}
    while step < steps:
        vf = max_vf*step/steps
        #print >>f, 0.0, " ", cryA.p(vf), " ", cryA.g(vf)
        manB_vf[step] = vf
        manB_p[step] = manB.p(vf)
        step += 1
    
    # Alder Crystal A
    cryA = AlderCrystal(1.0)
    f = open('bpd-crystal-a.dat','w')
    # Loop over the density:
    step = 1
    # Store EOS data:
    cryA_vf = {}
    cryA_p = {}
    print >>f, "# {alpha=",1.0,"} - vf, p, f, g, f_fluid, g_fluid"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, 0.0, " ", vf, " ", cryA.f(vf), " ", cryA.p(vf), " ", cryA.g(vf)
        cryA_vf[step] = vf
        cryA_p[step] = cryA.p(vf)
        step += 1
    f.close()
        
    # Alder Crystal B
    cryB = AlderCrystal(alpha)
    f = open('bpd-crystal-b.dat','w')
    # Loop over the density:
    max_vf = 0.7405
    step = 1
    # Store EOS data:
    cryB_vf = {}
    cryB_p = {}
    print >>f, "# {alpha=",alpha,"} - vf, p, f, g, f_fluid, g_fluid"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, 1.0, " ", vf, " ", cryB.f(vf), " ", cryB.p(vf), " ", cryB.g(vf)
        cryB_vf[step] = vf
        cryB_p[step] = cryB.p(vf)
        step += 1
    f.close()

    # Mansoori Fluid AB
    # Fixed composition line, NVT ensemble
    vf = 0.01
    max_vf = 0.7
    steps = 500
    f = open('bpd-Mansoori-fluid.dat','w')
    man = Mansoori(alpha,xB)
    while vf < max_vf:
        # Actually needs to be a pressure that gives 0.5*(vfa+vfb) = vf.
        pres = man.p(vf)
        while True:
            delta = 1e-8
            vfa = manA.vf(pres+delta)
            vfb = manB.vf(pres+delta)
            vftd = 0.5*(vfa+vfb)
            
            vfa = manA.vf(pres)
            vfb = manB.vf(pres)
            vft = 0.5*(vfa+vfb)
            #print vft," ",vftd
            if m.fabs(vft-vf) < 1.0e-6:
            #    print "break"
                break
            pres = pres + (vf-vft)/((vftd-vft)/delta)
            #vf = vf - (self.p(vf)-pres)/((self.p(vf+delta)-self.p(vf))/delta)
            # Phase-seperated free energy:
            fps = (xA*manA.f(vfa)+xB*manB.f(vfb))
            
        print >>f, vf, " ", man.Z(vf)," ", man.f(vf), " ", fps, " ", (man.f(vf)-fps), " ", (0.5*(vfa+vfb))-vf, " ", pres, " ", manA.p(vfa)," ",manB.p(vfb)
        vf += max_vf/steps
    f.close()
    
    # Mansoori AB Fluid, looking for phase seperation.
    nu1 = 0.5
    nu2 = 0.1
    nuR = 0
    while nuR <= 1:
        nu1c = nuR*nu1
        nu1f = nu1 - nu1c
        mann = Mansoori(alpha,nu2/(nu1f+nu2))
        u1f = mann.uA(nu1f+nu2)
        if nu1c > 0.4:
            u1c = cryA.u(nu1c)
        else:
            u1c = 0.0
        print nu1c, " ", u1f," ",u1c," "
        nuR += 1.0/200
    
    
# xB get TRASHED below this point ----------------------------------------------------


#    xB = 0.0
#    max_xB = 1.0
#    f = open('bpd-fluid-ab.dat','w')
#    fc = open('bpd-crystal-ab.dat','w')
#    max_vf = 0.8
#    print >>f, "# {alpha=",alpha,"} vf x_B PV/NkT P [s_A/kT] f [Nkt] g[NkT] S/Nk G/NkT"
#    while xB < max_xB+0.0000001:
#        man = Mansoori(alpha, xB)
#        # Loop over the density:
#        step = 1
#        while step <= steps:
#            vf = max_vf*step/steps
#            #if man.g(vf) < (1-xB)*cryA.g(cryA.vf(man.p(vf)))+ xB*cryB.g(cryB.vf(man.p(vf))):
#            pres = man.p(vf)
#            if pres > 11.930:
#              vfc = cryA.vf(pres)
#              Ag = cryA.g(vfc)
#              Af = cryA.f(vfc)
#            else:
#              vfc = manA.vf(pres)
#              Ag = manA.g(vfc)
#              Af = manA.f(vfc)
#              
#            if pres > 11.930/(alpha*alpha*alpha):
#              vfc = cryB.vf(pres)
#              Bg = cryB.g(vfc)
#              Bf = cryB.f(vfc)
#            else:
#              vfc = manB.vf(pres)
#              Bg = manB.g(vfc)
#              Bf = manB.f(vfc)
#            
#            #if man.g(vf) > Ag and man.g(vf) > 0.0:
#            #    print >>f, xB, " ", man.p(vf), " ", man.g(vf)
#            g = (1.0-xB)*Ag+xB*Bg
#            fh = (1.0-xB)*Af+xB*Bf
#            #if man.f(vf) > fh and man.f(vf) > 0.0:
#            print >>f, xB, " ", pres, " ", man.g(vf)
#            print >>fc, xB, " ", pres, " ", g
#            step += 1
#        print >>f, " "
#        print >>fc, " "
#        xB += 1.0/steps
#    f.close()
#    fc.close()

    # Examine a constant pressure line
    pres = 20.0
    prex_cx = 8.15
    #pres = pres_cx*sqrt(2.0)
    print "Switch pressures at ",(pres_cx*sqrt(2))," and ",(pres_cx*sqrt(2)/(alpha*alpha*alpha))
    f = open('bpd-p11.930.dat','w')
    fpres = open('bpd-p11.930-pres.dat','w')
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
    
    max_xB = 1.0
    steps = 500.0
    xB = max_xB/steps
    print >>f, "# {alpha=",alpha,"} vf x_B PV/NkT P [s_A/kT] f [Nkt] g[NkT] S/Nk G/NkT"
    while xB < max_xB:
        man = Mansoori(alpha, xB)
        vfab = man.vf(pres)
        manN = Mansoori(1.0/alpha, (1.-xB))
        #vfabN = manN.vf(pres/(alpha*alpha*alpha))
        #vfabN = manN.vf(pres)
        vfabN = vfab
        Ndel = 1.0
        N = 100000.0
        N2 = xB*N
        N1 = N-N2
        xBa = 1.0-((N1)/(N1+N2))
        xBb = 1.0-((N1+Ndel)/(N1+N2+Ndel))
        xBc = (N2+Ndel)/(N1+N2+Ndel)
        manA = Mansoori(alpha, xBa )
        manB = Mansoori(alpha, xBb )
        manC = Mansoori(alpha, xBc )
        vfabA = manA.vf(pres)
        vfabB = manB.vf(pres)
        vfabC = manC.vf(pres)
        chempotA = man.g(vfab) + (N1+N2+Ndel)*(manB.g(vfabB)-man.g(vfab))/(Ndel)
        chempotB = man.g(vfab) + (N1+N2+Ndel)*(manC.g(vfabC)-man.g(vfab))/(Ndel)
        #chempotA = manA.f(vfabA) + (N1+N2+Ndel)*(manB.f(vfabA)-manA.f(vfabA))/(Ndel)
        #chempotB = manA.f(vfabA) + (N1+N2+Ndel)*(manC.f(vfabA)-manA.f(vfabA))/(Ndel)
	    #
        print >>fpres, xB, " ", cryA.p(vfa), " ", manB.p(vfb), " ", man.p(vfab), " ", manN.p(vfabN), " ", alpha, " ", (1.-xB)," ", xB
#        print >>f, xB, " ", ((1-xB)*Ag+xB*Bg), " ", man.g(vfab), " ", man.uA(vfab), " ", man.uB(vfab), " ", cryA.g(vfa), " ", ((1-xB)*Af+xB*Bf), " ", man.f(vfab)," ", vfab, " ", vfa," ",vfb, cryA.p(vfa), " ", manB.p(vfb), " ", man.p(vfab)
#        print >>f, xB, " ", ((1-xB)*Ag+xB*Bg)," ", cryA.g(vfa), " ", cryA.f(vfa), " ", man.g(vfab)," ",  man.uA(vfab), " ", man.uB(vfab), " ", chempotA, " ", chempotB, " ", ((1.0-xB)*cryA.g(vfa)+xB*Bg)
        print >>f, xB, " ", cryA.g_abs(vfa), " ", man.g(vfab)," ", ((1.0-xB)*(cryA.g(vfa))+xB*Bg), " ", man.uA(vfab), " ", cryA.u(vfa)
#        print >>f, xB, " ", vfa," ", cryA.f(vfa), " ", vfab," ", man.f(vfab)
        xB += 1.0/steps
    f.close()
    fpres.close()
    
    # Pure A - Pure B ties lines (Seperate fluids, Crystal_A+Fluid_B, Crystal_A+Crystal_B)
    # I know Pa = 11.930, and so that Pb = 11.930/alpha^3
    # TODO This has moved a bit!  Pa now = pres_cx*sqrt(2)
    xB = 1.0 - 1.0/14.0
    man = Mansoori(alpha, xB)
    vf = 0.55
    conv = pi*((1.0-xB)+xB*alpha**3)/(6.*vf)
    print "Man AB conv=",conv
    print "Man AB   =",xB," vf=",vf," p=",man.p(vf), " g=", man.g(vf), " P\\rho/kT=",(man.p(vf)*conv)," Z=",man.Z(vf)
    vf = 0.55
    vfb = 0.315
    vfac = 0.563
    conv = (pi*((1.0-xB)+xB*alpha**3)/(6.*vf))
    print "Man AB conv=",conv
    print "Man AB   =",xB," vf=",vf," p=",man.p(vf), " g=", man.g(vf), " P\\rho/kT=",(man.p(vf)*conv)," Z=",man.Z(vf)
    print "Man B  vf=",vfb," p=",manB.p(vfb), " g=", manB.g(vfb)
    print "Ald A  vf=",vfac," p=",cryA.p(vfac)," g=",cryA.g(vfb)
    print "ManAld  G=",((1-xB)*cryA.g(vfb) + xB*manB.g(vfb))
    
    
    
    # Cry_A-to-Fluid_B Tie lines?
    xB = 0.0
    max_xB = 1.0
    f = open('bpd-tie-a-b.dat','w')
    max_vf = 0.70
    print >>f, "# {alpha=",alpha,"} vf x_B PV/NkT P [s_A/kT] f [Nkt] g[NkT] S/Nk G/NkT"
    # Loop over the density:
    step = 1
    while step <= steps:
        vf = 0.4+(max_vf-0.4)*step/steps
        print >>f, 0.0, " ", cryA.p(vf), " ", cryA.g(vf)
        vfb = manB.vf(cryA.p(vf))
        if( vfb < 0.495 ):
            print >>f, 1.0, " ", manB.p(vfb), " ", manB.g(vfb)
        else:
            vfc = cryB.vf(cryA.p(vf))
            print >>f, 1.0, " ", cryB.p(vfc), " ", cryB.g(vfc)
            
        print >>f, " "
        step += 1
    xB += 1.0/steps
    f.close()
    
    # Following Sollich, construct Helmholtz free energy surfaces in the partial-volumes plane
    max_nu = 0.5
    steps = 60.0
    f = open('bpd-Fpartial.dat','w')
    fcry = open('bpd-Fpartial-cry.dat','w')
    nu1 = 0.0
    while nu1 < max_nu:
        nu2 = 0.0
        while nu2 < max_nu:
            #  Construct the phase for this point:
            vf = nu1 + nu2
            if( nu1 == 0 and nu2 == 0 ):
                xB = 0.0
            else:
                xB = nu2/(nu1+nu2)
                
            # Look at the two crystals - full phase seperation
            cryA = AlderCrystal(1.0)
            cryB = AlderCrystal(alpha)
            fc = 0
            if nu1 > 0.52:
                fc += (1.0-xB)*cryA.f(nu1)
            else:
                fc += (1.0-xB)*manA.f(nu1)
            if nu2 > 0.52:
                fc += xB*cryA.f(nu2)
            else:
                fc += xB*manB.f(nu2)
            print >>fcry, nu1, " ", nu2, " ", fc
            
            # Build the fluid:
            man = Mansoori(alpha,xB)
            if( vf < 0.74 ):
                print >>f, nu1, " ", nu2, " ", fc-man.f(vf)
                
            nu2 += max_nu/steps
        print >>f, ""
        print >>fcry, ""
        nu1 += max_nu/steps
    f.close()
    fcry.close()
    
    # Open files for coex curves:
    fa = open('bpd-coex-a.dat','w')
    
    # Now, loop over pressures.
    max_pres = 50.0
    steps = 40
    step = 0
    pres = pres_cx*sqrt(2.0)
    mans = Mansoori(alpha, xB)
    while pres < max_pres:
        # In each case, N-R over the composition, looking for the cross-over for pure A etc.
        vf_cry = cryA.vf(pres)
        Gcry = cryA.g_abs(vf_cry)
        xB_t = mans.xB_for_coex_pres(pres,Gcry,0.0)
        if step == 0:
            xB_t = 0.0
        mans = Mansoori(alpha,xB_t)
        vf_mans = mans.vf(pres)
            
        print pres, "/", max_pres, " ", xB_t
        print >>fa, xB_t, " ", pres, " ", vf_cry, " ", vf_mans, " ", (1.0-xB_t)*vf_mans/((1.0-xB_t)+xB_t*alpha**3.0), " ", (xB_t*alpha**3.0)*vf_mans/((1.0-xB_t)+xB_t*alpha**3.0)
        fa.flush()
        # Loop increment:
        step = step + 1
        pres = pres + max_pres/(1.0*steps)
    
    # Close output files:
    fa.close();
        