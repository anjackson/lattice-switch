#!/usr/bin/env python
# --- --- ---
#
# Used to calculated useful factors for building
# lattices with the same number of sites, but different
# proportions.
#
# Arg can be the total number of sites, or the system dimensions
# as "nx ny nz".
#
# Will output all the common factors.  the last one should be the 
# most 'square'
#
# --- --- ---
# Andy Jackson - 18/Jul/2006
# $Id$
# --- --- ---

# Dependencies:
import sys

# Parse Command line arguments:
if len(sys.argv) == 2:
    n = int(sys.argv[1])
elif len(sys.argv) == 4:
    n = int(sys.argv[1]) * int(sys.argv[2]) * int(sys.argv[3])
else:
    print "usage: 3Dfactors.ny <N>"
    print "usage: 3Dfactors.ny <Nx> <Ny> <Nz>"
    sys.exit()
    

# Status:
print "Finding 3 integers that multiply up to %i" % n

# Hunt down all the factors [brute force method]:
i = 0
while i <= n:
    j = i
    while j < n:
        k = j
        while k < n:
            if (i*j*k) == n:
                print "Got %i %i %i = %i" %(i, j, k, n)
            k+=1
        j+=1
    i+=1

print "Done."
