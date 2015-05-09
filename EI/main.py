#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
############  Main script that coordinates  ##############
#############  all parts of the program   ################
##########################################################

from Scripts.rungen import *
from Scripts.io import *
import argparse
import sys
import os
import random
import shutil

############ GLOBALS DEFINITION ############
installdir="/home/timis/rdklol/EI/" # Installation directory
rundir = installdir+'Calc/'     # Jobs directory
###### END GLOBALS DEFINITION ##############

                
if __name__ == "__main__":
    # Parse commandline arguments
    ss = "\n***********************************************************"
    ss += "\n******** Welcome to SIMPLIFY! Let's get started. **********\n"
    ss += "***********************************************************"
    print ss
    # parse command line arguments
    parser = argparse.ArgumentParser()
    args = parsecommandline(parser,installdir)
    # make sure user specifies input
    if len(sys.argv) < 2 :
        print "\nYou didn't specify any arguments. Please try again. Exiting..\n"
        exit()
    # parse input file
    if (args.infile):
        parseinput(args)
    # check for jobs directory
    rundir = args.rundir+'/' if (args.rundir) else rundir
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    ################### START MAIN ####################    
    if (args.rgen):        
        ligs = readdict(installdir+'/Ligands/ligands.dict') # get ligands
        # unique single ligands
        lnums = random.sample(range(0, len(ligs)), min([len(ligs),int(args.rgen[0])]))
        for k in range(0,len(lnums)):
            args.lig = []
            args.ligocc = []
            args.lig.append(ligs.keys()[lnums[k]])
            rungen(installdir,rundir,args)
        dnum = int(args.rgen[0]) - len(ligs)
        if (dnum>0):
            lnums0 = random.sample(range(0, len(ligs)), dnum)
            lnums1 = random.sample(range(0, len(ligs)), dnum)
            for k in range(0,dnum):
                args.lig = []
                args.ligocc = []
                args.lig.append(ligs.keys()[lnums0[k]])
                args.lig.append(ligs.keys()[lnums1[k]])
                rungen(installdir,rundir,args) # run structure generation
    else:
        rungen(installdir,rundir,args) # run structure generation
    ss = "***********************************************************"
    ss += "\n***** Thank you for using SIMPLIFY. Have a nice day! ******\n"
    ss += "***********************************************************"
    print ss
    #################### END MAIN #####################
    
    

