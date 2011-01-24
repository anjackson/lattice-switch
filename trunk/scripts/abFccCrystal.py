#!/usr/bin/env python
#
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
    steps = 500
    max_vf = 0.7405
  
    # Alder Crystal A
    cryA = AlderCrystal(1.0, 0.0)
    cryB = AlderCrystal(alpha, 1.0)
    f = open('bpd-crystal-ab-fcc.dat','w')
    # Loop over the density:
    step = 1
    # Store EOS data:
    print >>f, "# {alpha=",1.0,"} - vf, p, f, g, f_fluid, g_fluid"
    while step < steps:
        vf = 0.5+(max_vf-0.5)*step/steps
        pres = cryA.p(vf)
        vf_b = cryB.vf(pres)
        vf_t = xA*vf + xB*vf_b
        f_t = xA*cryA.f(vf) + xB*cryB.f(vf_b)
        g_t = xA*cryA.g(vf) + xB*cryB.g(vf_b)
        if m.fabs(pres-cryB.p(vf_b)) < 1e-5:
            print >>f, vf_t, pres, f_t, cryA.f(vf)-f_t, g_t, cryA.g(vf)
        step += 1
    f.close()
    