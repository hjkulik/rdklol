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
from random import randint
import time
import pickle
import sys
import copy
import random
from math import sqrt

############ GLOBALS DEFINITION ############
installdir="/home/timis/rdklol/EI/" # Installation directory

# Temporary smiles dictionary of handwritten SMILES strings.
# Want to put database support for chembl/emolecules in here shortly.
smilesdict={'benzoate': 'c1ccc(cc1)C(=O)[O-]',
            'myristate': 'CCCCCCCCCCCCCC(=O)[O-]',
            'acetate': 'CC(=O)[O-]'}
            
ancores={'bicarbonate':'bicarbonate.mol','bisulfate':'bisulfate.mol',
         'bisulfite':'bisulfite.mol', 'dihydrogenphosphate':'dihydrogenphosphate.mol',
         'dihydrogenphosphite':'dihydrogenphosphite.mol','nitrate':'nitrate.mol',
         'oxalate':'oxalate.mol','perchlorate':'perchlorate.mol'}

licores={'ethyl':('ethyl.mol',1),'methyl':('methyl.mol',0),'carbonyl':('carbonyl.mol',0),'carbonyl':('carbonyl.mol',0),
            'cyanide':('cyanide.mol',0),'hydroxyl':('hydroxyl.mol',0),'amine':('amine.mol',0),
            'nitro':('nitro.mol',1),'pyridine':('pyridine.mol',3),'PPh3':('PPh3.mol',2),
            'benzene':('benzene.mol',4),'thiocyanate':('thiocyanate.mol',1),
            'acetonitrile':('acetonitrile.mol',1)} # ligands and connecting atoms indices

mcores={'ferrocene':('ferrocene_core.mol','PDR.mol')} # cores and characteristing groups
###### END GLOBALS DEFINITION ##############

def getligs():
    a=''
    for key in licores:
        a+=key+' '
    return a

def getanions():
    a=''
    for key in ancores:
        a+=key+' '
    return a
    
def getcores():
    a=''
    for key in mcores:
        a+=key+' '
    return a

def distance(R1,R2):
    dx = R1[0] - R2[0] 
    dy = R1[1] - R2[1] 
    dz = R1[2] - R2[2] 
    d = sqrt(dx**2+dy**2+dz**2)
    return d

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
    userinput = userinput.split('.mol')[0] 
    userinput = userinput.split('_core')[0]
    """ Check if user input matches defined dictionary. Could add support for new user defined cores. Not for now.
    """
    if userinput not in mcores:
        print("We didn't find the core structure: %s in the dictionary. Try again!\n" %(userinput))
        exit()
    if not glob.glob(installdir+'Cores/'+mcores[userinput][0]):
        print("We can't find the core structure file %s%s right now! Something is amiss.\n" %(installdir+'Cores/',mcores[userinput][0]))
        exit()
    else:
        # load core mol file (without hydrogens)
        core=Chem.MolFromMolFile(installdir+'Cores/'+mcores[userinput][0],removeHs=False)
        core.cg=mcores[userinput][1]
    return core
    
def lig_load(userinput):
    userinput = userinput.split('.mol')[0] 
    if userinput not in licores:
        print("We didn't find the ligand structure: %s in the dictionary. Try again!\n" %(userinput))
        exit()
    if not glob.glob(installdir+'Ligands/'+licores[userinput][0]):
        print("We can't find the ligand structure file %s%s right now! Something is amiss.\n" %(installdir+'Ligands/',licores[userinput][0]))
        exit()
    else:
        # load lig mol file (with hydrogens)
        lig=Chem.MolFromMolFile(installdir+'Ligands/'+licores[userinput][0],removeHs=False)
        lig.atID=licores[userinput][1]
    return lig

def anion_load(userinput):
    userinput = userinput.split('.mol')[0]+'.mol' 
    if not glob.glob(installdir+'Anions/'+userinput):
        print("We can't find the anion structure file %s%s right now! Something is amiss.\n" %(installdir+'Anions/',userinput))
        exit()
    else:
        # load anion mol file (without hydrogens)
        anion=Chem.MolFromMolFile(installdir+'Anions/'+userinput,removeHs=False)
    return anion    
    
def add_lig(core,lig,cgmol,cgmolH):
    # functionalizes core with ligands
    atom1 = 0   
    atoms =[]  
    maxH = 0
    hrem = 0
    Hrems = []
    # loop over atoms to find connection atom (most hydrogens connected)
    for i,atom in enumerate(cgmolH.GetAtoms()):
        numH = 0 
        for bond in atom.GetBonds():
            a = bond.GetBeginAtom()
            e = bond.GetEndAtom()
            if (a.GetIdx()!=i and a.GetAtomicNum()==1):
                numH+=1 
                hrem = a.GetIdx()
            elif(e.GetIdx()!=i and e.GetAtomicNum()==1):
                numH+=1
                hrem = e.GetIdx()
        if (numH > maxH):
            atoms = []
            Hrems = []
            atoms.append(i)
            maxH = numH
            Hrems.append(hrem)
        elif(numH==maxH):
            atoms.append(i)
            Hrems.append(hrem )
    # get random atom    
    idx = randint(0,len(atoms)-1)
    atom1 = atoms[idx]
    Hremove = Hrems[idx]
    # replace hydrogen by new ligand
    efm = Chem.EditableMol(cgmolH)    # get editable mol
    efm.RemoveAtom(Hremove)    
    fmol = efm.GetMol() # get normal mol back    
    atom2 = lig.atID+fmol.GetNumAtoms()   # connecting atom from ligand
    # combine cg with ligand
    mol = AllChem.CombineMols(fmol,lig)        
    eefm = Chem.EditableMol(mol)
    eefm.AddBond(atom1,atom2,Chem.BondType.SINGLE)            
    mol = eefm.GetMol() # get normal mol back
    Chem.SanitizeMol(mol)
    AllChem.EmbedMolecule(mol) # get coordinates    
    AllChem.MMFFOptimizeMolecule(mol) # optimize structure
    # combine with core structure
    ats0 = core.GetSubstructMatches(cgmolH) # to be removed from core structure    
    ats1 = core.GetSubstructMatches(cgmol) # to be removed from core structure        
    ats2 = mol.GetSubstructMatches(cgmol)  # to be aligned with core structure
    core3D = convert_to_mol3D(core)
    matom = core3D.GetAtom(0) # metal center coordinates
    # align with core structure and combine
    for at in ats1:
        # find best reflection
        AllChem.AlignMol(mol,core,atomMap=zip(ats2[0],at))    
        mol03D = convert_to_mol3D(mol)
        d0 = distance(mol03D.centermass(),matom.coords()) # get distance from center
        AllChem.AlignMol(mol,core,atomMap=zip(ats2[0],at),reflect=True)    
        mol13D = convert_to_mol3D(mol)
        d1 = distance(mol13D.centermass(),matom.coords()) # get distance from center
        # check for correct alignment compared to metal center
        if (d0 > d1):
            AllChem.AlignMol(mol,core,atomMap=zip(ats2[0],at),reflect=True)    
        core = AllChem.CombineMols(core,mol)        
    ratoms = sorted([v for atom in ats0 for v in atom])[::-1] # get inverse sorted list of atoms to remove
    ecore = Chem.EditableMol(core)
    # remove duplicate atoms
    for atom in ratoms:
        ecore.RemoveAtom(atom)
    core = ecore.GetMol()    
    return [core,fmol]
    

def structgen(args,rootdir,ligand,ligoc):
    strfiles = []
    # load molecule core
    core = core_load(args.core)   
    if (ligoc != 0):
        # load characteristic group
        cgmolH=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=False)    
        cgmol=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=True)    
        for i in range(0,int(ligoc)):
            # load ligand with H2
            lig = lig_load(ligand)
            # functionalize core
            [core,cgmolH] = add_lig(core,lig,cgmol,cgmolH)
    # convert molecule to mol3D
    core3D = convert_to_mol3D(core)
    # load anion, add hydrogens and convert to mol3D
    anion = anion_load(args.anion)
    #anionHs = prepmol(anion)
    AllChem.EmbedMolecule(anion)
    an3D = convert_to_mol3D(anion)
    # get minimum distance between core and anion
    mindist = core3D.size + an3D.size
    # assign reference point
    Rp = core3D.centermass()
    # generate multiple geometric arrangements
    Nogeom = int(args.anionsnum) # number of different combinations
    # Generate base case (separated structures)
    base3D = protate(copy.deepcopy(an3D),Rp,[20*mindist,0.0,0.0])
    mols = []
    maxdist = mindist+float(args.maxd) # Angstrom, distance of non-interaction    
    mindist = mindist-float(args.mind) # Angstrom, distance of non-interaction    
    fname = rootdir+'/'+args.core[0:3]+args.anion[0:2] 
    for i in range(0,Nogeom):        
        # generate random sequence of parameters for rotate()
        R = random.uniform(mindist,maxdist)
        phi = random.uniform(0.0,360.0)
        theta = random.uniform(0.0,180.0)
        thetax = random.uniform(0.0,360.0)
        thetay = random.uniform(0.0,360.0)
        thetaz = random.uniform(0.0,360.0)
        # translate
        tr3D = protate(copy.deepcopy(an3D), Rp,[R,theta,phi])
        newmol = cmrotate(tr3D,[thetax,thetay,thetaz])
        newmol.writemxyz(core3D,fname+str(i))
        strfiles.append(fname+str(i))
    print 'Generated ',Nogeom,' structures!\n'
    return strfiles
