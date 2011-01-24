#!/usr/bin/env python
#
# Class to deal with parsing lattice-switch simulation log files.
#
# Anj
# $Id$
#

# Base on the UserDict class:
from UserDict import UserDict
# Uses regular expressions for matching:
import re
# Can use command line arguments:
import sys

class LSLogParser(UserDict):
    "Parses log files produced by the LS codes and allows useful into to be extracted."
    
    def __init__(self, logfile="hsls.out"):
        "Initialised using the name of the log file, which defaults to 'hsls.out'."
        # Initialise the UserDict superclass:
        UserDict.__init__(self)
        # Store the logfile location:
        self["logfile"] = logfile
        # parse the log file:
        self.parseLogFile()
        
    def parseLogFile(self):
        "Parses the density out of the log file, and the error too."
        f = open(self["logfile"],'r')
        # Build the RegEx:
        self.reg = {};
        self.reg["density0_line"] = re.compile("^ CV lat0 dens .*$");
        self.reg["vf0_line"] = re.compile("^ CV lat0 vf .*$");
        self.reg["coa0_line"] = re.compile("^ CV lat0 c/a .*$");
        self.reg["density1_line"] = re.compile("^ CV lat1 dens .*$");
        self.reg["vf1_line"] = re.compile("^ CV lat1 vf .*$");
        self.reg["coa1_line"] = re.compile("^ CV lat1 c/a .*$");
        self.reg["accr_line"] = re.compile("^ AR .*$");
        self.reg["df_line"] = re.compile("^ F BTIP BA\(Df\) Df.*$")
# ANJ No-errors form:
#        self.reg["df_line"] = re.compile("^ F TIP-VS Df.*$")
        # Loop through the file, line by line:
        line = f.readline()
        while line:
            # Remove trailing space:
            line = line.rstrip()
            # See if this is a line about the density, store if so:
            for rek in self.reg.keys():
              if self.reg[rek].match(line):
                self[rek] = line;
            # Go on to the next line:
            line = f.readline()
        # Analyse the stored lines:
        self["density0"] = self.splitLine("density0_line",4)
        self["density0_error"] = self.splitLine("density0_line",6)
        self["density1"] = self.splitLine("density1_line",4)
        self["density1_error"] = self.splitLine("density1_line",6)
        self["vf0"] = self.splitLine("vf0_line",4)
        self["vf0_error"] = self.splitLine("vf0_line",6)
        self["vf1"] = self.splitLine("vf1_line",4)
        self["vf1_error"] = self.splitLine("vf1_line",6)
        self["coa0"] = self.splitLine("coa0_line",4)
        self["coa0_error"] = self.splitLine("coa0_line",6)        
        self["coa1"] = self.splitLine("coa1_line",4)
        self["coa1_error"] = self.splitLine("coa1_line",6)        
        self["dr"] = self.splitLine("accr_line",4)
        try:
            self["dvol"] = self.splitLine("accr_line",7)
        except:
            pass
        if "df_line" in self:
            # F BTIP BA(Df) Df/N= -2.009336644725901 +/- 5.920946318591509E-4
            self["df"] = self.splitLine("df_line", 5)
            self["df_err"] = self.splitLine("df_line", 7)
# ANJ No-errors form:
#            self["df"] = self.splitLine("df_line", 4)
#            self["df_err"] = 0.0

    def splitLine(self, lkey, item):    
        if self.has_key(lkey):
            line = self[lkey]
            # Split on spaces and store the density and the error:
            bits = line.split(" ")
            return bits[item]
        else:
            return ""

    def getTwoKey(self, key):
        if self[key+"0"] and self[key+"1"]:
            return "["+self[key+"0"] + "," + self[key+"1"]+"]"
        elif self[key+"1"]:
            return self[key+"1"]            
        else:
            return self[key+"0"]
        
    def getTwoKeyErr(self, key):
        if self[key+"0_error"] and self[key+"1_error"]:
            return "["+self[key+"0_error"] + "," + self[key+"1_error"]+"]"
        elif self[key+"1_error"]:
            return self[key+"1_error"]            
        else:
            return self[key+"0_error"]
        
   
    def getDensity(self):
        return self.getTwoKey("density")
    
    def getDensityError(self):
        return self.getTwoKeyErr("density")
        
    def getVolFrac(self):
        return self.getTwoKey("vf")
    
    def getVolFracError(self):
        return self.getTwoKeyErr("vf")
        
    def getCOverA(self):
        if self["coa1"]:
            return self["coa1"]
        else:
            return self["coa0"]
    
    def getCOverAError(self):
        if self["coa1"]:
            return self["coa1_error"]
        else:
            return self["coa0_error"]
        
    def getDR(self):
        return self["dr"]
        
    def getDVol(self):
        return self["dvol"]
    
    def getDf(self):
        return self["df"]
        
    def getDfErr(self):
        return self["df_err"]
        
# TODO Should check for WARNING, ERROR, MELTED, ETC and refuse if there were problems.
    

# If this script is invoked as a stand-alone program:
if __name__ == '__main__':
    # Use the first command line argument as the log file name:
    lslp = LSLogParser(sys.argv[1])
    # Attempt to parse the density:
    print "Density: "+lslp.getDensity() + " +/- " + lslp.getDensityError()
    print "Volume Fraction: "+lslp.getVolFrac() + " +/- " + lslp.getVolFracError()
    print "c/a: "+lslp.getCOverA() + " +/- " + lslp.getCOverAError()
    print "dr: "+lslp.getDR()
    print "dVol: "+lslp.getDVol()
    
            