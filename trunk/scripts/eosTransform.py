#!/usr/bin/env python
#
# Combines data with Mansoori and Alder to produce co-ex curve data.
#

import sys
from math import log,pi,fabs
from AlderCrystal import *
    
# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Choose the file to load:
    if len(sys.argv) >= 2:
        dff = sys.argv[1]
    else:
        print "usage: eosTransform.py <eos.dat> [<alpha>] [<xB>]";
        sys.exit()
        
    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) >= 3:
        alpha = float(sys.argv[2])
    else:
        alpha = 1.0
        
    # Choose the fraction that is B:
    if len(sys.argv) >= 4:
        nB = float(sys.argv[3])
        if nB < 0.0:
            xB = - nB
        else:
            xB = nB/(1+nB)
        max_xB = xB        
    else:
        # Scan the whole range of xB
        xB = 0.0
        max_xB = 1.0
        
    # Echo
    print '# alpha = ', alpha, ' xB = ', xB

    # Set up the reference crystal.        
    cryA = AlderCrystal(alpha)
    
    # Calculate the close-packing VF (cheating for now)
    # FCC
    if nB == 0:
        vf_cp = 0.7404804895

    # AB(CsCl) \alpha = 0.73
    if nB == 1:
        vf_cp = 0.7272876005

    # AB2 \alpha = 0.58
    if nB == 2:
        vf_cp = 0.7763822606

    # AB13 \alpha = 0.58
    if nB == 13:
        vf_cp = 0.7128773506
#        vf_cp = 0.73
    
    print '# vf_cp = ',vf_cp

    # Open the file and loop over the lines
    dfp = open(dff,'r')
    # Loop through the file, line by line:
    line = dfp.readline()
    while line:
        # Remove space from both ends:
        line = line.strip()
        # Pres dens StdErr(dens) vf StdErr(vf) c/a StdErr(c/a) dR dVol
        # Skip comments:
        if line[0] != '#':
            bits = line.split()
            e_p = float(bits[0])
            e_d = float(bits[1])
            e_derr = float(bits[2])
            e_vf = float(bits[3])
            e_vferr = float(bits[4])
            e_coa = float(bits[5])
            e_coaerr = float(bits[6])
            e_dr = float(bits[7])
            e_dv = float(bits[8])
            
            # Volume:
            vol = (pi/(6.0*e_vf)) * ((1.0-xB)+xB*(alpha**3.0))
            
            # Compressability - PV/NkT
            Z = e_p*vol
            
            # Reduced volume - 1
            rvol = e_vf/vf_cp
            rvol = (1-rvol)/rvol
            
            # Further reduce:
            #Z = Z - 3/rvol  
            
            print rvol," ", Z, " ", e_vferr

        # Next line:
        line = dfp.readline()
