#!/usr/bin/env python
#
# Tester for the various equations drawn from the literature.
#
# The aim it to test the codes by reproducing the relevant figures shown in the original papers.
#

import sys
import math as m
from CSFluid import *
from AlderCrystal import *
from Mansoori import *

class VolFracLooper:
    "Loop over the volume fraction and output the given parameters"
    pass

        
# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) == 2:
        alpha = float(sys.argv[1])
    else:
        alpha = 1.0
        
    csf = CSFluid(alpha)
    # Loop over the density:
    steps = 500
    max_vf = 0.7405
    f = open('cs-eos.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - vf, Z, p"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, vf, " ", csf.Z(vf), " ", csf.p(vf)
        step += 1
    f.close()
        
    # Generate pV_0/NkT v. V/V_0 graph from Alder paper.
    f = open('Alder-eos-V0.dat', 'w')
    alder = AlderCrystal(alpha)
    step = 1
    print >>f, "# {alpha=",alpha,"} - vf, Z, p"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, alder.Vred(vf), " ", alder.p(vf)/m.sqrt(2), " ", csf.p(vf)/m.sqrt(2)
        step += 1
    print >>f," "
    print >>f,"1.352 8.426 8.426"
    print >>f,"1.495 8.426 8.426"
    f.close()
    
    # Generate S(wrtIG) v. V_0/V table from Alder paper.  Also allow comparision with Mansoori
    man = Mansoori(1.0,0.0)
    f = open('Alder-SIG-V0.dat', 'w')
    print >>f, "# {alpha=",alpha,"} - V_0/V S vf"
    dens = 0.0001
    densS = 0.01
    densM = 0.65
    while dens <= densM:
        vf = dens*pi*sqrt(2.0)/6.0
        print >>f, dens, " ", -alder.f_fluid(vf), " ", -csf.f(vf), " ", -man.f(vf),  " ", man.g(vf)
        dens = dens + densS
    f.close()
    
    # Check the one data point from Alder:
    Vr = 1.35
    vf = pi*sqrt(2.0)/(6.0*Vr)
    print " Check Vr=",Vr," vf=",vf," S-Sig=-6.01 == ",-alder.f(vf)
    
    # Generate Helmholtz v. VolFrac graph (Alder Crys and CS Fluid)
    f = open('Alder-AvVf.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - vf, F_c, F_f"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, alder.Vred(vf), alder.f(vf), " ", alder.f_fluid(vf), " ", csf.f(vf)
        step += 1
    f.close()
        
    # Generate Gibbs v. Pressure graph
    f = open('Alder-GvP-CrysFluid.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - P, G_c, vf"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, alder.p(vf), " ", alder.g(vf), " ", alder.p_fluid(vf), " ", alder.g_fluid(vf), " ", vf
        step += 1
    f.close()
    
    # And the fluid:
    f = open('CS-GvP-Fluid.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - P, G_f, vf"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, csf.p(vf), " ", csf.g(vf), " ", vf
        step += 1
    f.close()
    
    # Generate the Gibbs v. Pressure graph for Mansoori v. FCC_A and FCC_B
    man = Mansoori(1.0,0.0)
    f = open('Mansoori-FCCa-GvP.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - P, G_f"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, man.p(vf), " ", man.g(vf), " ", csf.g(vf)
        step += 1
    alder = AlderCrystal(1.0)
    step = 1
    print >>f, ""
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, alder.p(vf), " ", alder.gM(vf)
        step += 1
    f.close()
    # PRESSURE AND CHEMICAL POTENTIAL C.F. Bartlett [1990]
    man = Mansoori(1.0,0.0)
    alder = AlderCrystal(1.0)
    f = open('Mansoori-FCCa-UvVf.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - P, G_f"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, vf, " ", man.p(vf), " ", alder.p(vf), " ", man.uA(vf), " ", alder.u(vf)
        step += 1
    vf = 0.4954
    print "At Fluid vf=",vf,", P=",man.p(vf)/sqrt(2)," uA=",man.uA(vf)
    print "At Fluid vf=",vf,", f=",man.f(vf)," g=",man.g(vf)
    vf = 0.54776
    print "At Crystal vf=",vf,", P=",alder.p(vf)/sqrt(2)," u=",alder.u(vf)
    print "At Crystal vf=",vf,", f=",alder.f(vf)," g=",alder.g(vf)
    f.close()
    # Alpha = 0.58
    man = Mansoori(0.58,1.0)
    f = open('Mansoori-FCCb-GvP.dat', 'w')
    step = 1
    print >>f, "# {alpha=",alpha,"} - P, G_f"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, man.p(vf), " ", man.g(vf)
        step += 1
    alder = AlderCrystal(0.58)
    step = 1
    print >>f, ""
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print >>f, alder.p(vf), " ", alder.gM(vf)
        step += 1
    f.close()

    # Generate S(wrtIG) v. V_0/V table from Alder paper.  Also allow comparision with Mansoori
    alpha = 1.0/3.0
    xB = 0.5
    man = Mansoori(alpha,xB)
    f = open('Mansoori-d3-x0.5.dat', 'w')
    print >>f, "# {alpha=",alpha," x_B=",xB,"} "
    dens = 0.0001
    densS = 0.01
    densM = 0.65
    while dens <= densM:
        vf = dens*pi*sqrt(2.0)/6.0
        print >>f, vf, " ", man.f(vf), " ", csf.f(vf), " ", alder.f_fluid(vf), " ", man.g(vf), " ", man.GH(vf)
        dens = dens + densS
    f.close()
    
    alpha = 1.0
    xB = 0.0
    man = Mansoori(alpha,xB)
    alder = AlderCrystal(alpha)
    f = open('Mansoori-d1-x0.0.dat', 'w')
    print >>f, "# {alpha=",alpha," x_B=",xB,"} "
    dens = 0.0001
    densS = 0.01
    densM = 0.65
    while dens <= densM:
        vf = dens*pi*sqrt(2.0)/6.0
        print >>f, vf, " ", man.f(vf), " ", csf.f(vf), " ", alder.f_fluid(vf), " ", man.g(vf), " ", man.GH(vf)
        dens = dens + densS
    f.close()
    print "Density of Fluid 0.49546, want P=11.930: ", csf.p(0.49546), " ", man.p(0.49546)
    print "Density of Crystal 0.54784, want P=11.930: ", alder.p(0.54784)
