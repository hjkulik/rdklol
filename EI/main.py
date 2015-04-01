#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
##########  Top level script that coordinates  ###########
##########    generation of structues, input   ###########
##########         files, jobscripts           ###########
##########################################################

from Scripts.structgen import structgen
from Scripts.tcgen import tcgen
from Scripts.jobgen import jobgen
import argparse
import sys
import os
import shutil

if __name__ == "__main__":
    # Parse commandline arguments
    ss = "\n***********************************************************"
    ss += "\n********* Welcome to rdklol! Let's get started. ***********\n"
    ss += "***********************************************************"
    print ss
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--core", help="core structure type") #e.g. ferrocene
    parser.add_argument("-a","--anion", help="anion type") #e.g. bisulfate, nitrate, perchlorate -> For binding
    parser.add_argument("-na","--anionsnum", help="number of anion copies") #different geometric arrangements for calculating binding energy
    parser.add_argument("-lc","--ligandcore", help="ligand core structure type") #e.g. bidentate, monodentate -> Functionalize ferrocene
    parser.add_argument("-l","--ligand", help="ligand structure name") #e.g. acetate (in smilesdict) -> Functionalize ferrocene
    parser.add_argument("-smi","--smiles", help="(optional) manual SMILES string") #e.g 'O=C[O-]'
    parser.add_argument("-q","--dbquery", help="(optional) database to query") #database to query
    parser.add_argument("-n","--dbnum", help="(optional) database list-based number to query", type=int) #list element in database
    parser.add_argument("-rn","--regnum", help="(optional) database reg number to query", type=int) #reg number of ligand in database (ligkey)
    parser.add_argument("-f","--force", help="(optional) force build of structure with bad non-bonded distances", action="store_true")
    args=parser.parse_args()
    # make sure user specifies input
    if len(sys.argv) < 2 :
        print "\nYou didn't specify any arguments. Please try again. Exiting..\n"
        exit()
    #### START MAIN ####
    # create folder for runs
    rootdir = 'Calcs/'+args.core+'_'+args.anion
    if os.path.isdir(rootdir):
        flagdir=raw_input('\nDirectory '+rootdir +' already exists. Remove? (y/n) :')
        if (flagdir=='y' or flagdir=='Y'):
            shutil.rmtree(rootdir)
            os.mkdir(rootdir)
    else:
        os.mkdir(rootdir)
    # generate structures
    structgen(args,rootdir)
    # generate terachem input files
    
    # generate jobfiles
    
    
    
    
    ss = "***********************************************************"
    ss += "\n****** Thank you for using rdklol. Have a nice day! *******\n"
    ss += "***********************************************************"
    print ss
    
    

