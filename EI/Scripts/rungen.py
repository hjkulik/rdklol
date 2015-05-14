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
import itertools

def randomgen(installdir,rundir,args):
    ligs = readdict(installdir+'/Ligands/ligands.dict') # get ligands
    # load global variables
    licores = readdict(installdir+'/Ligands/ligands.dict')
    # get all combinations of ligands
    combos = []
    for i in range(1,7):
        for combs in itertools.combinations(range(0,len(ligs)),i):
            combos.append(combs)
    # get a sample of these combinations
    samps = random.sample(range(0,len(combos)),int(args.rgen[0]))
    # loop over samples
    for c in samps:
        combo = combos[c]
        args.lig = []
        args.ligocc = []
        args.totocc = 0 
        combol = list(combo)
        random.shuffle(combol)
        print combol
        for cj in combol:
            rocc = random.randint(1,6-args.totocc+1)
            trocc = rocc*int(len(licores[ligs.keys()[cj]][1:]))
            if args.totocc+trocc <= 6:
                args.lig.append(ligs.keys()[cj])
                args.ligocc.append(rocc)
                args.totocc += trocc
        if len(args.lig) > 0 :
                rungen(installdir,rundir,args) # run structure generation
    return args

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
    skip = False
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
            print "\nSkipping folder "+rootdir+" ...\n\n"
            skip = True
    else:
        os.mkdir(rootdir)
    if not skip:
        print 'Generating runs in folder '+rootdir+' ..'
        # generate structures
        strfiles = structgen(installdir,args,rootdir,ligands,ligocc)
        # generate terachem input files
        jobdirs = tcgen(args,strfiles,lig)
        print 'Created the terachem input files!'
        # generate jobfiles
        jobgen(args,jobdirs,lig)
        print 'Created the jobscripts!\n\n\n'
    
    

