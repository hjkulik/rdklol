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
import argparse
import re
import sys
import os
import shutil

def parseinput(args):
    for line in open(args.infile):
        li = line.strip()
        if not li.startswith("#") and len(li)>0: # remove comments/empty lines
            l = line.split('#')[0] # remove comments
            l = filter(None,re.split(' |,|\t',l))
            if (l[0]=='-c'):
                args.core = l[1]
            if (l[0]=='-a'):
                args.anion = l[1]
            if (l[0]=='-na'):
                args.anionsnum = l[1]
            if (l[0]=='-m'):
                args.method = l[1]
            if (l[0]=='-b'):
                args.basis = l[1]
            if (l[0]=='-sp'):
                args.spin = l[1]
            if (l[0]=='-ch'):
                args.charge = l[1]
            if (l[0]=='-q'):
                args.queue = l[1]
            if (l[0]=='-g'):
                args.gpus = l[1]
            if (l[0]=='-maxd'):
                args.maxd = l[1]
            if (l[0]=='-mind'):
                args.mind = l[1]
            if (l[0]=='-lig'):
                args.lig = l[1:]
            if (l[0]=='-suff'):
                args.suff = l[1].strip('\n')
            if (l[0]=='-dispersion'):
                args.dispersion = l[1].strip('\n')
                
if __name__ == "__main__":
    # Parse commandline arguments
    ss = "\n***********************************************************"
    ss += "\n******** Welcome to SIMPLIFY! Let's get started. **********\n"
    ss += "***********************************************************"
    print ss
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile",help="specified in input file")
    parser.add_argument("-c","--core", help="core structure with currently available: "+getcores()) #e.g. ferrocene
    parser.add_argument("-a","--anion", help="anion type with currently available: "+getanions()) #e.g. bisulfate, nitrate, perchlorate -> For binding
    parser.add_argument("-na","--anionsnum", help="number of anion copies") #different geometric arrangements for calculating binding energy
    parser.add_argument("-lig","--lig", help="ligand structure name with currently available: "+getligs()) #e.g. acetate (in smilesdict) -> Functionalize ferrocene
    parser.add_argument("-smi","--smiles", help="(optional) manual SMILES string") #e.g 'O=C[O-]'
    parser.add_argument("-dbq","--dbquery", help="(optional) database to query") #database to query
    parser.add_argument("-n","--dbnum", help="(optional) database list-based number to query", type=int) #list element in database
    parser.add_argument("-rn","--regnum", help="(optional) database reg number to query", type=int) #reg number of ligand in database (ligkey)
    parser.add_argument("-f","--force", help="(optional) force build of structure with bad non-bonded distances", action="store_true")
    parser.add_argument("-maxd","--maxd", help="Maximum distance above cluster size for molecules placement maxdist=size1+size2+maxd", action="store_true")
    parser.add_argument("-mind","--mind", help="Minimum distance above cluster size for molecules placement mindist=size1+size2+mind", action="store_true")
    # terachem arguments    
    parser.add_argument("-xyz","--xyzfile", help="Input file")
    parser.add_argument("-m","--method", help="electronic structure approach for terachem job. Specify UHF/UDFT to enforce levelshifting and unrestricted calculation for OS singlets (default: b3lyp).")
    parser.add_argument("-b","--basis", help="basis for terachem job (default: LACVP*).") 
    parser.add_argument("-s","--spin", help="spin for system (default: singlet).")
    parser.add_argument("-ch","--charge", help="charge for system (default: neutral).") 
    parser.add_argument("-x","--extra", nargs="+", help="extra arguments in syntax tckeyword=value. suggestions: nstep (opt), maxit (scf), scf=diis+a, min_tolerance>4.5e-4/min_tolerance_e>1.0e-6, min_coordinates=cartesian, orbitalswrtfrq=<num>,multibasis=<basis file>.")
    parser.add_argument("-tc","--terachem", help="full path to custom terachem installation.")
    # jobscript arguments
    parser.add_argument("-mem","--memory", help="memory reserved per thread for job file (default: 2G.")
    parser.add_argument("-t","--time", help="time requested for queueing system (default: 168hrs).")
    parser.add_argument("-q","--queue", help="queue wildcards (default: [dg]*).")
    parser.add_argument("-g","--gpus", help="number of GPUS (default: 1).")    
    parser.add_argument("-j","--job", help="prefix for pdb/xyz structure file and jobname.")
    parser.add_argument("-dispersion","--dispersion", help="dispersion yes or no.")
    parser.add_argument("-suff","--suff", help="suffix for jobs folder.")
    args=parser.parse_args()
    # make sure user specifies input
    if len(sys.argv) < 2 :
        print "\nYou didn't specify any arguments. Please try again. Exiting..\n"
        exit()
    if (args.infile):
        parseinput(args)
    #### START MAIN ####
    # create folder for runs    
    numlig=-1
    if (args.lig):
        ligands = []
        ligocc=[]
        for entry in args.lig:
            if not entry.isdigit():
                ligands.append(entry) # get ligands
                ligocc.append('1') # initialize ligands occurence
                numlig += 1
            else:
                ligocc[numlig] = entry # get occurrence in the molecule
    else:
        ligands=='0'
        ligocc=0
    for i,lig in enumerate(ligands):
        lig = False if ligands=='0' else lig
        occur = ligocc[i]
        if (args.suff):
            rootdir = 'Calcs/'+args.core[0:4]+lig+occur+'_'+args.anion+args.suff
        else:
            rootdir = 'Calcs/'+args.core[0:4]+lig+occur+'_'+args.anion
        if os.path.isdir(rootdir):
            flagdir=raw_input('\nDirectory '+rootdir +' already exists. Remove? (y/n) : ')
            if (flagdir=='y' or flagdir=='Y' or flagdir=='Yes' or flagdir=='yes'):
                shutil.rmtree(rootdir)
                os.mkdir(rootdir)
        else:
            os.mkdir(rootdir)
        print 'Generating runs in folder '+rootdir+' ..\n'
        # generate structures
        strfiles = structgen(args,rootdir,lig,occur)
        # generate terachem input files
        jobdirs = tcgen(args,strfiles,lig)
        print 'Created the terachem input files!\n'
        # generate jobfiles
        jobgen(args,jobdirs,lig)
        print 'Created the jobscripts!\n'
    ss = "***********************************************************"
    ss += "\n***** Thank you for using SIMPLIFY. Have a nice day! ******\n"
    ss += "***********************************************************"
    print ss
    
    

