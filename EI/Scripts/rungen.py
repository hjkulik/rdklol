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
from collections import Counter
import itertools

def randomgen(installdir,rundir,args):
    # load global variables
    if (args.core == 'ferrocene'):
        licores = readdict(installdir+'/Ligands/ligands_ferro.dict')
    else:
        licores = readdict(installdir+'/Ligands/ligands.dict')
    # get all combinations of ligands
    combos = []
    for i in range(1,7):
        for combs in itertools.combinations(range(0,len(licores)),i):
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
        for cj in combol:
            rocc = random.randint(1,6-args.totocc+1)
            trocc = rocc*int(len(licores[licores.keys()[cj]][1:]))
            if args.totocc+trocc <= 6:
                args.lig.append(licores.keys()[cj])
                args.ligocc.append(rocc)
                args.totocc += trocc
        if len(args.lig) > 0 :
                rungen(installdir,rundir,args) # run structure generation
    return args

def counterSubset(list1, list2):
        c1, c2 = Counter(list1), Counter(list2)
        for k, n in c1.items():
            if n > c2[k]:
                return False
        return True

def getconstsample(no_rgen,args,licores):
    samp = []
    # 3 types of constraints: ligand, ligocc, coord
    # get ligand and ligocc
    get = False
    occup=[]
    coord = int(args.coord) if (args.coord) else 6
    combos = []
    generated = 0 
    for i in range(1,coord+1):
        combos += (list(itertools.combinations_with_replacement(range(0,len(licores)),i)))
    if (args.lig):
        licindx = []
        for idl,lig in enumerate(args.lig):
            locc = args.ligocc[idl] if (args.ligocc and len(args.ligocc) > idl) else 1
            occup.append(locc)
            # get index in licores
            for ii,li in enumerate(licores):
                if (lig in li):
                    for a in range(0,int(occup[idl])):
                        licindx.append(ii)
        # get right combos
        random.shuffle(combos)
        for combo in combos:
            if (counterSubset(licindx,combo)):
                # get total denticity
                totdent = 0
                for l in combo:
                    totdent += int(len(licores[licores.keys()[l]][2:]))
                if totdent == coord:
                    # add combo
                    samp.append(combo)
                    generated += 1
            if (generated >= no_rgen):
                break
    else:
        for combo in combos:
            # get total denticity
            totdent = 0
            for l in combo:
                totdent += int(len(licores[licores.keys()[l]][2:]))
            if totdent == coord:
                # add combo
                samp.append(combo)
                generated += 1
            if (generated >= no_rgen):
                break
    return samp

    
def constrgen(installdir,rundir,args):
    # load global variables
    if (args.core == 'ferrocene'):
        licores = readdict(installdir+'/Ligands/ligands_ferro.dict')
    else:
        licores = readdict(installdir+'/Ligands/ligands.dict')
    # get a sample of these combinations
    samps = getconstsample(int(args.rgen[0]),args,licores)
    # loop over samples
    for combo in samps:
        args.lig = []
        args.ligocc = []
        for cj in set(combo):
            lcount = Counter(combo)
            rocc = lcount[cj]
            args.lig.append(licores.keys()[cj])
            args.ligocc.append(rocc)
        rungen(installdir,rundir,args) # run structure generation
    return args

def rungen(installdir,rundir,args):
    # check for specified ligands/functionalization
    ligocc = []
    if (args.lig):
        ligands = args.lig
        if (args.ligocc):
            ligocc = args.ligocc
        else:
            ligocc.append('1')
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
        strfiles,charge = structgen(installdir,args,rootdir,ligands,ligocc)
        # generate terachem input files
        jobdirs = tcgen(args,strfiles,lig,charge)
        print 'Created the terachem input files!'
        # generate jobfiles
        jobgen(args,jobdirs,lig)
        print 'Created the jobscripts!\n\n'
    
    

