#!/usr/bin/env python
#
# Basic numerical implementation of the Carnahan-Starling equation for a pure HS fluid.
#
# REF: The Journal of Chemical Physics -- Volume 51, No. 2, pp.635-636
#

import sys
from math import log,sqrt,pi,e

class CSFluid:
    """
    The Carnahan-Starling Equation.
    ---
    Length unit is \sigma_A.
    """
    
    def __init__(self, alpha=1.0):
        "Initialise parameters for the single crystal."
        # Store the crystal unit-scaling parameters:
        if alpha == 1.0:
            # If alpha is 1.0, then we are a pure A crystal:
            self.xB = 0.0
            self.alpha = 1.0
        else:
            # Otherwise, this is a pure B crystal with reduced radius alpha:
            self.xB = 1.0
            self.alpha = alpha
                
    def V(self, vf):
        "Calculate the volume [1/N.sigma_A] for this volume fraction."
        return (pi*((1.0-self.xB)+self.xB*(self.alpha**3.0)))/(6.0*vf)
    
    def Vred(self,vf):
        "Calculate the reduced volume used in the original paper."
        return sqrt(2.0)*self.V(vf)
    
    def Z(self,vf):
        "Return the PV/NkT according to the Carhahan-Startling equation for the hard-sphere fluid."
        return (1.0+vf+vf**2.0-vf**3.0)/((1.0-vf)**3.0)

    def p(self, vf):
        "Returns the pressure at the given volume fraction."
        return self.Z(vf) / self.V(vf)
    
    def f(self,vf):
        "Returns the free-energy wrt IG at the given volume fraction."
        return (3.0-2.0*vf)/((1.0-vf)**2.0) - 3.0
    
    def g(self,vf):
        "Returns the Gibbs free energy wrt IG at pressure for volfrac."
        return self.f(vf) + self.Z(vf) - 1.0 - log(self.Z(vf))

    def g_abs(self,vf):
        "Return the Gibbs free energy (absolute), in units of [NkT]."
        return self.f(vf) + self.Z(vf) - 1.0 + log(self.V(vf)) 
        
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
    step = 1
    print "# {alpha=",alpha,"} - vf, p, f, g, f_fluid, g_fluid"
    while step < steps:
        vf = 0.3+(max_vf-0.3)*step/steps
        print vf, " ", csf.Z(vf), " ", csf.p(vf)
        step += 1
