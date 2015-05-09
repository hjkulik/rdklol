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

from Scripts.structgen import *
from Scripts.tcgen import tcgen
from Scripts.jobgen import jobgen
from Scripts.io import *
import argparse
import sys
import os
import shutil

def rungen(installdir,rundir,args):
    # check for specified ligands/functionalization
    if (args.lig):
        ligands = args.lig
        ligocc = args.ligocc
        for i in range(len(ligocc),len(ligands)):
            ligocc.append('1')
        lig = ''.join("%s%s" % (l[0:3],ligocc[i]) for i,l in enumerate(ligands)) # folder name
    else:
        ligands =[]
        lig = ''
        ligocc = ''
    # create folder for runs and check if it already exists
    if (args.anion):
        rootdir = rundir+args.core[0:4]+lig+args.anion[0:4]
    else:
        rootdir = rundir+args.core[0:4]+lig
    if (args.suff):
        rootdir += args.suff
    if os.path.isdir(rootdir):
        flagdir=raw_input('\nDirectory '+rootdir +' already exists. Replace? (y/n) : ')
        if (flagdir=='y' or flagdir=='Y' or flagdir=='Yes' or flagdir=='yes'):
            shutil.rmtree(rootdir)
            os.mkdir(rootdir)
        else:
            print "\nProgram exiting..."
            exit(0)
    else:
        os.mkdir(rootdir)
    print 'Generating runs in folder '+rootdir+' ..\n'
    # generate structures
    strfiles = structgen(installdir,args,rootdir,ligands,ligocc)
    # generate terachem input files
    jobdirs = tcgen(args,strfiles,lig)
    print 'Created the terachem input files!\n'
    # generate jobfiles
    jobgen(args,jobdirs,lig)
    print 'Created the jobscripts!\n'
    
    

