#!/usr/bin/env python
#
# This python script runs NPT simulations at a range of pressures,
# and determines the density and c/a ratio for each pressure.
#
# $Id$
# ----
# Andrew N Jackson, 29/08/2006.
#

import string
import math as m
import sys

# Use my log parser:
from LSLogParser import LSLogParser

# Reference crystal: Alder FCC:
from AlderCrystal import AlderCrystal
from Mansoori import Mansoori

# Class to handle building parameter filed:
class ParFile:
    "This class builds simulation parameter files from templates."

    def __init__(self, template):
        " Initialise the template."
        self.templ = open(template).read()
        # Initialise an empty dictionary:
        self.subs = dict()
        # Default substitutions:
        self.subs['volfrac'] = 0.6
        self.subs['l0cOa'] = 1.0
        self.subs['l1cOa'] = 1.0
        self.subs['dr'] = 0.01
        self.subs['t_tot'] = 100000
        self.subs['t_equib'] = 5000
        self.subs['t_out'] = 1000
        self.subs['t_check'] = 2000
        self.subs['ss_new_tpm'] = 'y'
        self.subs['ss_read_weights'] = 'n'
        self.subs['evolve_weights'] = 'n'
        self.subs['reset_u'] = 'n'
        
    def setPar(self, key, value):
        "Sets the paremeter 'key' to the supplied value."
        self.subs[key] = value
        
        
    def createParameterFile(self, template='hsls.p%{pressure)s.df.par'):
        "Puts the parameters into a file."
        # Create the parameter file from the template:
        pars = self.templ % self.subs
        # Create a filename (based on another template:
        parf = open((template % self.subs),'w')
        parf.write(pars)
        parf.close()
        
    def setForWFEvol(self):
        self.subs['t_tot']   = 6000000
        self.subs['t_check'] =   10000
        self.subs['evolve_weights'] = 'y'
        self.subs['reset_u'] = 'y'
    
    def setForDfRun(self):
        self.subs['t_tot']   = 10000000
        self.subs['t_check'] =  1000000
        self.subs['evolve_weights'] = 'n'
        self.subs['reset_u'] = 'n'


# ---------------------------------------------------------------

# Base on the UserDict class:
from UserDict import UserDict

# Class to read EOS files:
class EOSFile(UserDict):
    "Reads equation of state files in as a dict indexed on pressure."
    
    def __init__(self, eosfile="eos.dat"):
        # Initialise the UserDict superclass:
        UserDict.__init__(self)
        # Read in from the file:
        fin = open(eosfile)
        for line in fin:
            # Skip comments:
            if line.lstrip().startswith('#'): continue
            # Split the line into its conponents:
            pres, dens, dens_e, vf, vf_e, cOa, cOa_e, dR, dVol = line.split()
            par = dict()
            par['pres'] = pres
            par['dens'] = dens
            par['vf'] = vf
            par['cOa'] = cOa
            par['dR'] = dR
            par['dVol'] = dVol
            #Store it:
            self[int(pres)] = par
    
# ---------------------------------------------------------------

import os

# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Default to NOT a dummy run:
    dummy_run = 0

    # Choose the alpha from CLARGS, no arg implies alpha = 1.0 (pure A crystal).
    if len(sys.argv) == 3 or len(sys.argv) == 4:
        alpha = float(sys.argv[1])
        nB = float(sys.argv[2])
        xB = nB/(1+nB) 
        print "#Got alpha=",alpha," x_B=",xB
    else:
        print "I need two arguments! <Alpha> <#Bs for each A> [-dummy-run]"
        sys.exit()
    # Pick up the dummy-run flag if it's there:
    if len(sys.argv) == 4:
        # If there is a 4th argument, this is a dummy run that should only analyse.
        dummy_run = 1
            
    # The classpath to use:
    #classpath = '../../bin';
    classpath = '/Home/ajackso1/Projects/LatticeSwitch/bin';
    #classpath = '/home/andrewj/Projects/LatticeSwitch/build';
    
    # The EOS file to use:
    eos_file = "eos.dat"
    
    # Command prefix - set to 'echo ' for dummy runs:
    if dummy_run == 1:
        cmd_pre = "echo "
    else:
        cmd_pre = ""
    
    # Print out results header:
    print "# Pres vf c/a f_ex g_ex df stderr(df)"
    
    # The template to use:
    nvt_template = './hsls.nvt.par.template'
    # Create a parameter file creator:
    parm = ParFile(nvt_template)
    
    # Open the EOS file:
    eos = EOSFile(eos_file)
    keys = eos.keys(); keys.sort()
    for pres in keys:
        #print " I Creating parameters for P=%(pres)d" % { 'pres':pres } 
        parm.setPar('pressure', pres)
        parm.setPar('volfrac', eos[pres]['vf'])
        parm.setPar('l0cOa', eos[pres]['cOa'])
        parm.setPar('l1cOa', -1)
        parm.setPar('dr', eos[pres]['dR'])
        # Build the file id string:
        prestr = 'p%(pressure)d' % parm.subs
        presdir = "Df."+prestr
        # Create a directory and move into it:
        if not os.path.isdir(presdir): os.mkdir(presdir)
        os.chdir(presdir)
        #os.system("echo ' I PWD=$PWD'")
        # If no TPM file exists,do the w/f evol run:
        tpm = 'TPM.'+prestr
        wfoutfile = "hsls.p%(pr)s.wf.par.out" % { 'pr':pres }
        if not os.path.exists(wfoutfile):
            # Set up for weight-function run:
            parm.setForWFEvol()
            parm.createParameterFile('hsls.p%(pressure)s.wf.par')
            # Run this simulation:
            cmd = "java -Xmx1024m -cp %(cp)s net.anjackson.physics.ls.LSSim hsls.p%(pr)s.wf.par > %(wfo)s 2>&1" % { 'cp':classpath, 'pr':pres, 'wfo':wfoutfile }
            os.system(cmd_pre+cmd)
            # Copy TPM into TPM.in and into a file for storage:
            os.system(cmd_pre+"cp TPM.dat "+tpm)
            # Copy other files; pdf.*.dat, final.xyz
            os.system(cmd_pre+"cp pdf.vs.dat pdf.vs.wf.dat")
            os.system(cmd_pre+"cp pdf.tpm.dat pdf.tpm.wf.dat")
            os.system(cmd_pre+"cp final.wrl final.wf.wrl")
        
        # If there is a TPM file, but no results file for this run:
        outfile = 'hsls.'+prestr+'.df.par.out'
        if os.path.exists(tpm):
            if not os.path.exists(outfile):
                # Copy TPM file into place:
                os.system(cmd_pre+"cp "+tpm+" TPM.in")
                # Set up for dF-eval run:
                parm.setForDfRun()
        
                # create a parameter using the given filename template:    
                parm.createParameterFile('hsls.p%(pressure)s.df.par')
        
                # Run this simulation:
                cmd = "java -Xmx1024m -cp %(cp)s net.anjackson.physics.ls.LSSim hsls.p%(pr)s.df.par > %(out)s 2>&1" % { 'cp':classpath, 'pr':pres, 'out':outfile }
                os.system(cmd_pre+cmd)

        # If both result files exist:
        if os.path.exists(tpm) and os.path.exists(outfile):
            # Parse the output file
            lslp = LSLogParser(outfile)
        
            # Set up the reference crystal
            cry = AlderCrystal()

            # Print results: p dens vf c/a df
            vf =  float(eos[pres]['vf'].strip())
            vol = (m.pi/(6.0*vf))*((1.0-xB)+xB*(alpha**3))
            print pres," ", vf, " ", eos[pres]['cOa']," ", (cry.fM(vf)+float(lslp.getDf()))," ", ((cry.fM(vf)+float(lslp.getDf()))+pres*vol-1.0-m.log(pres*vol)), " ", lslp.getDf(), " ", lslp.getDfErr(), " ", ((cry.fM(vf)+float(lslp.getDf()))+pres*vol-1.0+m.log(vol))

        # Change back to the parent directory:
        os.chdir('..')
