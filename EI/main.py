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
from Classes.globalvars import *
from Classes.mol3D import mol3D
from Classes.atom3D import atom3D
import argparse
from math import sqrt
from math import floor
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
        randomgen(installdir,rundir,args)
    else:
        rungen(installdir,rundir,args) # run structure generation
    ss = "***********************************************************"
    ss += "\n***** Thank you for using SIMPLIFY. Have a nice day! ******\n"
    ss += "***********************************************************"
    print ss
    #################### END MAIN #####################
    
    

