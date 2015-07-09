#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
########## This script handles input/output  #############
##########################################################

# import std modules
from rdkit import Chem
import glob
import os
import re
import argparse
import sys

def parseinput(args):
    args.ox = 0
    for line in open(args.infile):
        li = line.strip()
        if not li.startswith("#") and len(li)>0: # remove comments/empty lines
            l = line.split('#')[0] # remove comments
            l = filter(None,re.split(' |,|\t',l))
            if (l[0]=='-c'):
                args.core = l[1].lower()
            if (l[0]=='-a'):
                args.anion = l[1].lower()
            if (l[0]=='-na'):
                args.anionsnum = l[1]
            if (l[0]=='-m'):
                args.method = l[1].lower()
            if (l[0]=='-b'):
                args.basis = l[1]
            if (l[0]=='-sp'):
                args.spin = l[1]
            if (l[0]=='-ch'):
                args.charge = l[1]
            if (l[0]=='-q'):
                args.queue = l[1].lower()
            if (l[0]=='-g'):
                args.gpus = l[1]
            if (l[0]=='-maxd'):
                args.maxd = l[1]
            if (l[0]=='-mind'):
                args.mind = l[1]
            if (l[0]=='-coord'):
                args.coord = l[1]
            if (l[0]=='-lig'):
                args.lig = [ll.lower() for ll in l[1:]]
            if (l[0]=='-ligocc'):
                args.ligocc = l[1:]
            if (l[0]=='-dispersion'):
                args.dispersion = l[1].strip('\n').lower()
            if (l[0]=='-rgen'):
                args.rgen = l[1:]
            if (l[0]=='-suff'):
                args.suff = l[1].strip('\n')
            if (l[0]=='-rdir'):
                args.rundir = l[1].strip('\n')
                if (args.rundir[-1]=='/'):
                    args.rundir = args.rundir[:-1]
            if (l[0]=='-ff'):
                args.ff = l[1].lower()
            if (l[0]=='-ox'):
                args.ox = int(l[1])
            if (l[0]=='-post'):
                args.postp = True
                

def parsecommandline(parser,installdir):
    parser.add_argument("-i","--infile",help="specified in input file")
    parser.add_argument("-c","--core", help="core structure with currently available: "+getcores(installdir)) #e.g. ferrocene
    parser.add_argument("-a","--anion", help="anion type with currently available: "+getanions(installdir)) #e.g. bisulfate, nitrate, perchlorate -> For binding
    parser.add_argument("-na","--anionsnum", help="number of anion copies") #different geometric arrangements for calculating binding energy
    parser.add_argument("-lig","--lig", help="ligand structure name with currently available: "+getligs(installdir)) #e.g. acetate (in smilesdict) -> Functionalize ferrocene
    parser.add_argument("-ligocc","--ligocc", help="number of corresponding ligands") # e.g. 1,2,1
    parser.add_argument("-coord","--coord", help="Coordination such as 4,5,6") # coordination e.g. 6 (octahedral)
    parser.add_argument("-rgen","--rgen", help="number of random generated molecules, overwrites lig and ligcorr")
    parser.add_argument("-smi","--smiles", help="(optional) manual SMILES string") #e.g 'O=C[O-]'
    parser.add_argument("-dbq","--dbquery", help="(optional) database to query") #database to query
    parser.add_argument("-n","--dbnum", help="(optional) database list-based number to query", type=int) #list element in database
    parser.add_argument("-rn","--regnum", help="(optional) database reg number to query", type=int) #reg number of ligand in database (ligkey)
    parser.add_argument("-f","--force", help="(optional) force build of structure with bad non-bonded distances", action="store_true")
    parser.add_argument("-maxd","--maxd", help="Maximum distance above cluster size for molecules placement maxdist=size1+size2+maxd", action="store_true")
    parser.add_argument("-mind","--mind", help="Minimum distance above cluster size for molecules placement mindist=size1+size2+mind", action="store_true")
    parser.add_argument("-rdir","--rundir",help="Directory for jobs",action="store_true")
    parser.add_argument("-ff","--ff",help="Force field optimize ligands",action="store_true")
    parser.add_argument("-postp","--postp",help="Post process results",action="store_true")
    # terachem arguments    
    parser.add_argument("-xyz","--xyzfile", help="Input file")
    parser.add_argument("-m","--method", help="electronic structure approach for terachem job. Specify UHF/UDFT to enforce levelshifting and unrestricted calculation for OS singlets (default: b3lyp).")
    parser.add_argument("-b","--basis", help="basis for terachem job (default: LACVP*).") 
    parser.add_argument("-s","--spin", help="spin for system (default: singlet).")
    parser.add_argument("-ch","--charge", help="charge for system (default: neutral).") 
    parser.add_argument("-ox","--ox", help="oxidation state for metal (default: neutral).") 
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
    return args
    
def readdict(fname):
    d = dict()
    f = open(fname,'r')
    txt = f.read()
    lines = filter(None,txt.splitlines())
    f.close()
    for line in lines:
        key = filter(None,line.split(':')[0])
        val = filter(None,line.split(':')[1])
        d[key] = filter(None,re.split(',| ',val))
    return d 
    
def getligs(installdir):
    licores = readdict(installdir+'/Ligands/ligands.dict')
    a=''
    for key in licores:
        a+=key+' '
    return a

def getanions(installdir):
    ancores = readdict(installdir+'/Anions/anions.dict')
    a=''
    for key in ancores:
        a+=key+' '
    return a
    
def getcores(installdir):
    mcores = readdict(installdir+'/Cores/cores.dict')
    a=''
    for key in mcores:
        a+=key+' '
    return a
    
def loaddata(fname):
    d = dict()
    f = open(fname)
    txt = f.read()
    lines = filter(None,txt.splitlines())
    for line in lines[1:]:
        l = filter(None,line.split(None))
        d[(l[0],l[1])] = [float(ll) for ll in l[2:]]
    f.close()
    return d

def loadcoord(installdir,coord):
    f = open(installdir+'/Data/'+coord+'.dat')
    txt = filter(None,f.read().splitlines())
    f.close()
    b = []
    for line in txt:
        l = filter(None,line.split(None))
        b.append([float(l[0]),float(l[1]),float(l[2])])
    return b
    

def core_load(installdir,userinput,mcores):
    userinput = userinput.split('.mol')[0] 
    userinput = userinput.split('_core')[0]
    """ Check if user input matches defined dictionary. Could add support for new user defined cores. Not for now.
    """
    if userinput not in mcores:
        print("We didn't find the core structure: %s in the dictionary. Try again!\nAvailable cores are:%s\n" %(userinput,getcores(installdir)))
        exit()
    if not glob.glob(installdir+'Cores/'+mcores[userinput][0]):
        print("We can't find the core structure file %s%s right now! Something is amiss.\n" %(installdir+'Cores/',mcores[userinput][0]))
        exit()
    else:
        # load core mol file (without hydrogens)
        core=Chem.MolFromMolFile(installdir+'Cores/'+mcores[userinput][0],removeHs=False)
        core.cg=mcores[userinput][0]
        core.cat = [int(l) for l in filter(None,mcores[userinput][1:])]
    return core
    
def lig_load(installdir,userinput,licores):
    userinput = userinput.split('.mol')[0] 
    userinput = userinput.split('.xyz')[0] 
    if userinput not in licores:
        print("We didn't find the ligand structure: %s in the dictionary. Try again!\nAvailable ligands are:%s\n" %(userinput,getligs(installdir)))
        exit()
    if not glob.glob(installdir+'Ligands/'+licores[userinput][0]):
        print("We can't find the ligand structure file %s%s right now! Something is amiss.\n" %(installdir+'Ligands/',licores[userinput][0]))
        exit()
    else:
        # load lig mol file (with hydrogens)
        flig = installdir+'Ligands/'+licores[userinput][0]
        lig = Chem.Mol()
        if ('.xyz' in flig):
            lig.f=flig
            lig.mol = False
        elif ('.mol' in flig):
            lig=Chem.MolFromMolFile(flig,removeHs=False)
            lig.mol = True
        lig.cat = [int(l) for l in licores[userinput][2:]]
        lig.charge = int(licores[userinput][1])
        lig.denticity = len(licores[userinput][2:])
    return lig

def anion_load(installdir,userinput,ancores):
    userinput = userinput.split('.mol')[0]+'.mol' 
    if not glob.glob(installdir+'Anions/'+userinput):
        print("We can't find the anion structure file %s%s right now! Something is amiss.\n" %(installdir+'Anions/',userinput))
        exit()
    else:
        # load anion mol file (without hydrogens)
        anion=Chem.MolFromMolFile(installdir+'Anions/'+userinput,removeHs=False)
        anion.charge = -1
    return anion   
    
    
     
    


