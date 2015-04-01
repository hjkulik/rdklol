#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script generates a collection  #############
########  of randomly placed anions around   #############
########   a functionalized ferrocene core   #############
######## and then creates the required input #############
######## and job files for running terachem  #############
########  calculations on these structures   #############
##########################################################

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalForceFields as FF
from rdkit.Chem import rdMolTransforms
# import custom modules
from tcgen import mybash
from Classes.mol3D import mol3D
from Classes.atom3D import atom3D
from rotate import protate
from rotate import cmrotate
# import std modules
import argparse
import psycopg2
import urllib2
import glob
import os
import shutil
import time
import pickle
import sys
import copy
import random

############ GLOBALS DEFINITION ############
installdir="/home/timis/rdklol/EI/" # Installation directory

# Temporary smiles dictionary of handwritten SMILES strings.
# Want to put database support for chembl/emolecules in here shortly.
smilesdict={'benzoate': 'c1ccc(cc1)C(=O)[O-]',
            'myristate': 'CCCCCCCCCCCCCC(=O)[O-]',
            'acetate': 'CC(=O)[O-]'}
            
licores={'bicarbonate':'bicarbonate.pdb','bisulfate':'bisulfate.pdb',
         'bisulfite':'bisulfite.pdb', 'dihydrogenphosphate':'dihydrogenphosphate.pdb',
         'dihydrogenphosphite':'dihydrogenphosphite.pdb','nitrate':'nitrate.pdb',
         'oxalate':'oxalate.pdb','perchlorate':'perchlorate.pdb'}

mcores={'ferrocene':('ferrocene_core.pdb',((1,2,3,4,5),(6,7,8,9,10)))}
###### END GLOBALS DEFINITION ##############


def convert_to_mol3D(rdmol):
    # create 3D molecule
    m3D = mol3D() 
    for i,atom in enumerate(rdmol.GetAtoms()):
        # get coordinates
        pos = rdmol.GetConformer().GetAtomPosition(i)
        # add atom to molecule
        m3D.addatom(atom3D(atom.GetSymbol(),[pos[0],pos[1],pos[2]]))
    return m3D

def core_load(userinput):
    userinput = userinput.split('.pdb')[0] 
    """ Check if user input matches defined dictionary. Could add support for new user defined cores. Not for now.
    """
    if userinput not in mcores:
        print("We didn't find the core structure: %s in the dictionary. Try again!\n" %(userinput))
        exit()
    if not glob.glob(installdir+'Cores/'+mcores[userinput][0]):
        print("We can't find the core structure file %s%s right now! Something is amiss.\n" %(installdir+'Cores/',mcores[userinput][0]))
        exit()
    else:
        # load core pdb file (without hydrogens)
        core=Chem.MolFromPDBFile(installdir+'Cores/'+mcores[userinput][0])
    return core

def anion_load(userinput):
    userinput = userinput.split('.pdb')[0]+'.pdb' 
    if not glob.glob(installdir+'Anions/'+userinput):
        print("We can't find the anion structure file %s%s right now! Something is amiss.\n" %(installdir+'Anions/',userinput))
        exit()
    else:
        # load anion pdb file (without hydrogens)
        anion=Chem.MolFromPDBFile(installdir+'Anions/'+userinput)
    return anion

def prepmol(mol):
    """Adds hydrogens to SMILES, turns it into a 3D"""
    molh=AllChem.AddHs(mol)
    AllChem.EmbedMolecule(molh)
    return molh

def structgen(args,rootdir):
    # check existence ferrocene molecule core
    core = core_load(args.core)   
    # Initialize core structure model
    coreHs = prepmol(core) # add H2 + convert to 3D 
    # convert molecule to mol3D
    core3D = convert_to_mol3D(coreHs)
    # load anion, add hydrogens and convert to mol3D
    anion = anion_load(args.anion)
    anionHs = prepmol(anion)
    an3D = convert_to_mol3D(anionHs)
    # get minimum distance between core and anion
    mindist = core3D.size + an3D.size
    # assign reference point
    Rp = core3D.centermass()
    # generate multiple geometric arrangements
    Nogeom = int(args.anionsnum) # number of different combinations
    # Generate base case (separated structures)
    base3D = protate(copy.deepcopy(an3D),Rp,[10*mindist,0.0,0.0])
    mols = []
    maxdist = 10 # Angstrom, distance of non-interaction    
    fname = rootdir+'/'+args.core+'_'+args.anion 
    for i in range(0,Nogeom):        
        # generate random sequence of parameters for rotate()
        R = random.uniform(0.8*mindist,mindist+maxdist)
        phi = random.uniform(0.0,360.0)
        theta = random.uniform(0.0,180.0)
        thetax = random.uniform(0.0,360.0)
        thetay = random.uniform(0.0,360.0)
        thetaz = random.uniform(0.0,360.0)
        # translate
        tr3D = protate(copy.deepcopy(an3D), Rp,[R,theta,phi])
        newmol = cmrotate(tr3D,[thetax,thetay,thetaz])
        newmol.writemxyz(core3D,fname+'_'+str(i))
    print 'Generated ',Nogeom,' structures!\n'
