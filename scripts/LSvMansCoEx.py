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
        print "usage: LSvMansCoEx.py <df.dat> [<alpha>] [<xB>]";
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
    print '# xB = ',xB

    # Set up the reference crystal.        
    cryA = AlderCrystal(alpha)
    man = Mansoori(alpha, xB)
    # Open the file and loop over the lines
    dfp = open(dff,'r')
    # Loop through the file, line by line:
    line = dfp.readline()
    # For interpolation:
    last_gabs = 0.0;
    last_p = -1.0;
    interp_max = 21.0;
    while line:
        # Remove space from both ends:
        line = line.strip()
        # Skip comments:
        if line[0] != '#':
            bits = line.split()
            df_p = float(bits[0])
            df_vf = float(bits[1])
            df_cOa = float(bits[2])
            df_f = float(bits[3])
            df_g = float(bits[4])
            df_df = float(bits[5])
            df_dferr = float(bits[6])
            df_gabs = float(bits[7])
            # Interpolate:
            interp = 0
            if( last_p < 0.0 ):
                interp = interp_max - 1
            while interp < interp_max:
                # Determined density of reference system.
                interp_factor = (interp_max - (interp+1.0))/interp_max
                pres = df_p - (df_p-last_p)*interp_factor
                vf = man.vf(pres)
                vf_cry = cryA.vf(pres)
                # Convert EX to ABS free energies...
#            Gcry = ((cryA.g_abs(vf_cry) - cryA.g(vf_cry)) + df_g) 
                Gcry = df_gabs - (df_gabs-last_gabs)*interp_factor 
                # For AB2
#                Gcry = Gcry + 2.0 - 0.15
                # For AB13
                Gcry = Gcry + 3.0 - 0.05
                if pres > 12:
                    print "# (G: ",pres, " " ,Gcry, " : ", cryA.g_abs(vf_cry), cryA.g(vf_cry), df_df
                    Xr_t = man.xB_for_coex_pres(pres, Gcry, nB)
                    xB_c = 1.0 - (xB/Xr_t)
                    xB_f = 1.0 - xB_c
                
                if pres <= 10:# or Xr_t < xB:
                    Xr_t = xB
                    xB_c = 0.0
                    xB_f = 1.0
            
                man_t = Mansoori(alpha, Xr_t)            
                vf_t = man_t.vf(pres)
            
                print vf, " ", xB, " ", pres, " ", Gcry, " ", man.g(vf), " ", man.g_abs(vf), " ", man_t.g_abs(vf), " ", (1.0-xB)*Gcry, " ", Xr_t, " ", xB_c, " ", (xB_c*Gcry + xB_f*man_t.g_abs(vf_t)), " ", vf_t*(1.0-Xr_t)/((1.0-Xr_t)+Xr_t*alpha**3.0), " ",vf_t*(Xr_t*alpha**3.0)/((1.0-Xr_t)+Xr_t*alpha**3.0)," ",df_vf*(1.0-xB)/((1.0-xB)+xB*alpha**3.0)," ",df_vf*(xB*alpha**3.0)/((1.0-xB)+xB*alpha**3.0)
                
                # Next interp point:
                interp = interp + 1
                
            # Store last line contents
            last_p = df_p
            last_gabs = df_gabs
            
        # Go on to the next line:
        line = dfp.readline()

    sys.exit()
    print "# {@alpha=",alpha,"} vf x_B PV/NkT P [s_A/kT] f [Nkt] g[NkT] g_abs [NkT]"
