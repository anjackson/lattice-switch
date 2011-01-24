#!/usr/bin/env python
#
# This python script runs NPT simulations at a range of pressures,
# and determines the density and c/a ratio for each pressure.
#
# $Id$
# ----
# Andrew N Jackson, 29/08/2006.
#

# Use my log parser:
from LSLogParser import LSLogParser 

# Class to handle building parameter filed:
class NPTParFile:
    "This class builds NPT simulation parameter files from templates."

    def __init__(self, template):
        " Initialise the template."
        self.templ = open(template).read()
        # Set up the basic parameter list:
        self.subs = {
    # initial volume fraction:
        'volfrac': '0.6',
    # Particle moves:
        'dr': '0.01',
        'particle.move.type': 'rw',
    # Pressure
        'pressure': '20',
        'dV': '0.001'
        }

    # Set the pressure:
    def setPressure(self,pres):
        self.subs['pressure'] = pres
        
    # function to write out parameter dictionary to a parameter file:
    def writePars(self, file):
        # Create the parameter file from the template:
        pars = self.templ % self.subs
        pfh = open( file, 'w')
        pfh.write(pars)
        pfh.close()

# Imports
import os

# ------------- Basic invocation parameters --------------

# Lomond classpath:
#classpath = '/home/andrewj/Projects/LatticeSwitch/build'
# Cluster classpath:
classpath = '/Home/ajackso1/Projects/LatticeSwitch/bin'

# Simulation Factory Invoker Class:
lsclass = 'net.anjackson.physics.ls.LSSim'

# Simulation parameter template file:
template = './hsls.npt.par.template'

# --------------------------------------------------------

# Init the parfile handler:
nptpar = NPTParFile(template)

# Output header:
print "# Pres dens StdErr(dens) vf StdErr(vf) c/a StdErr(c/a) dR dVol"

# Loop pressure values:
press = [ '200', '150', '120', '100', '90', '80', '70', '60', '50', '40', '35', '30', '28', '26', '24', '23', '22', '21', '20', '19', '18', '17', '16', '14', '12', '11', '10', '9', '8', '7', '6', '5', '4' ]
#press = [ '100', '90', '80', '70', '60', '50', '40', '35', '30', '28', '26', '24', '23', '22', '21', '20', '19', '18', '17', '16', '14', '12', '11', '10', '9', '8', '7', '6', '5', '4' ]

# Loop over the pressures:
for p in press:
    nptpar.setPressure(p)
    # Check if the output file already exists, continue if not:
    pfile = "hsls-p"+p+".par"
    ofile = pfile+".out"
    if not os.path.exists(ofile):
        nptpar.writePars(pfile)
        lssim = "java -cp "+classpath+" "+lsclass+" "+pfile+" > "+ofile+" 2>&1"
        os.system(lssim)
        # Copy the final WRL file to a unique name:
        wrlfile = "final-p"+p+".wrl"
        os.system("cp final.wrl "+wrlfile)
    
    # Extract the results from the output file.
    lslp = LSLogParser(ofile)
    print p + " " + lslp.getDensity() + " " + lslp.getDensityError() + \
       " " + lslp.getVolFrac() + " " + lslp.getVolFracError() + \
       " " + lslp.getCOverA() + " " + lslp.getCOverAError() + \
       " " + lslp.getDR() + " " + lslp.getDVol() 
