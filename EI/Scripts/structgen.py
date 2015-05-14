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
from geometry import *
from io import *
from Classes.globalvars import *
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
import numpy as np
from numpy import linalg as la

def distance(R1,R2):
    dx = R1[0] - R2[0] 
    dy = R1[1] - R2[1] 
    dz = R1[2] - R2[2] 
    d = sqrt(dx**2+dy**2+dz**2)
    return d

def setdiff(a,b):
    b = set(b)
    return [aa for aa in a if aa not in b]

def vecdiff(r1,r2):
    dr = [a-b for a,b in zip(r1,r2)]
    return dr

def convert_to_mol3D(rdmol):
    # create 3D molecule
    m3D = mol3D() 
    for i,atom in enumerate(rdmol.GetAtoms()):
        # get coordinates
        pos = rdmol.GetConformer().GetAtomPosition(i)
        # add atom to molecule
        m3D.addatom(atom3D(atom.GetSymbol(),[pos[0],pos[1],pos[2]]))
    return m3D

def rotation_params(r0,r1,r2):
    r10 = [a-b for a,b in zip(r1,r0)]
    r21 = [a-b for a,b in zip(r2,r1)]
    # angle between r10 and r32
    theta = 180*np.arccos(np.dot(r21,r10)/(la.norm(r21)*la.norm(r10)))/np.pi
    # get normal vector to plane r0 r1 r2
    u = np.cross(r21,r10)
    # check for collinear case
    if la.norm(u) < 1e-16:
        # pick random perpendicular vector
        if (abs(r21[0]) > 1e-16):
            u = [(-r21[1]-r21[2])/r21[0],1,1]
        elif (abs(r21[1]) > 1e-16):
            u = [1,(-r21[0]-r21[2])/r21[1],1]
        elif (abs(r21[2]) > 1e-16):
            u = [1,1,(-r21[0]-r21[1])/r21[2]]
    return theta,u

def add_lig(core,catom,lig,cgmol,cgmolH):
    #################################################
    ####### functionalizes core with ligands ########
    ############### for general cores ###############
    #################################################
    # get atom for connection
    numH = 0 
    hrem = 0 
    atom = cgmolH.GetAtomWithIdx(catom)
    for bond in atom.GetBonds():
        a = bond.GetBeginAtom()
        e = bond.GetEndAtom()
        if (a.GetAtomicNum()==1):
            numH += 1 
            hrem = a.GetIdx()
        elif (e.GetAtomicNum()==1):
            numH += 1 
            hrem = e.GetIdx()
    if (numH == 0):
        print "WARNING: Maximum valency for all atoms in molecule. Can't functionalize molecule further.."
        fmol = cgmolH
    else:
        # replace hydrogen by new ligand
        efm = Chem.EditableMol(cgmolH)    # get editable mol
        efm.RemoveAtom(hrem)    
        fmol = efm.GetMol() # get normal mol back    
        atom2 = lig.cat[0]+fmol.GetNumAtoms()   # connecting atom from ligand
        # combine cg with ligand
        mol = AllChem.CombineMols(fmol,lig)        
        eefm = Chem.EditableMol(mol)
        eefm.AddBond(catom,atom2,Chem.BondType.SINGLE)            
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
            AllChem.AlignMol(mol,core,atomMap=zip(ats2[0],at),reflect=True)    
            mol13D = convert_to_mol3D(mol)
            d0 = mol03D.mindist(core3D)
            d1 = mol13D.mindist(core3D)
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
    
    
def mcomplex(core,ligs,ligoc,MLbonds,installdir,licores):
    #################################################
    ####### functionalizes core with ligands ########
    ############## for metal complexes ##############
    #################################################
    coordbasef=['oct','tbp','thd','tri','li','one']
    metal = core.GetAtomWithIdx(0).GetSymbol()
    occs0 = []
    dentl = []
    toccs = 0 
    # find correct occurence for each ligand
    for i,lig in enumerate(ligs):
        dent_i = int(len(licores[lig][1:]))
        oc_i = int(ligoc[i]) if i < len(ligoc) else 1
        occs0.append(0)
        dentl.append(dent_i)
        for j in range(0,oc_i):
            if (toccs+dent_i <= 6):
                occs0[i] += 1
            toccs += dent_i
    # sort by descending denticity (needed for adjacent connection atoms)
    indcs = [i[0] for i in sorted(enumerate(dentl), key=lambda x:x[1],reverse=True)]    
    ligands = [ligs[i] for i in indcs]
    occs = [occs0[i] for i in indcs]
    coord = min(toccs,6) # complex coordination
    if toccs > 6:
        print "WARNING: coordination greater than 6, octahedral complex will be generated"
    # load base core for coordination
    corexyz = loadcoord(installdir,coordbasef[6-coord])
    # create molecule and add metal and base
    m3D = mol3D() 
    m3D.addatom(atom3D(metal,corexyz[0])) ## add metal
    core3D = mol3D() 
    core3D.addatom(atom3D(metal,corexyz[0])) ## add metal
    mcoords = core3D.GetAtom(0).coords()
    for m in range(1,coord+1):
        m3D.addatom(atom3D('X',corexyz[m])) ## add termination atoms
    # loop over ligands
    totlig = 0 
    for i,ligand in enumerate(ligands):
        for j in range(0,occs[i]):
            denticity = len(licores[ligand][1:])
            if not(ligand=='x' or ligand =='X') and (totlig-1+denticity < coord):
                lig = lig_load(installdir,ligand,licores) # load ligand
                lig3D = convert_to_mol3D(lig) # get 3D ligand
                catoms = lig.cat # connecting atom
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0
                ####################################################
                ## optimize geometry by minimizing steric effects ##
                ####################################################      
                if (denticity == 1):
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(m3D,lig3D.GetAtom(atom0),m3D.GetAtom(totlig+1))
                    r1 = lig3D.GetAtom(atom0).coords()
                    r2 = lig3D.centermass() # center of mass
                    rrot = r1
                    theta,u = rotation_params(r0,r1,r2) 
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,rrot,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,rrot,u,theta-180)
                    d2 = distance(mcoords,lig3D.centermass())
                    d1 = distance(mcoords,lig3Db.centermass())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # rotate around axis of symmetry and get best orientation
                    if lig3D.natoms > 2 :
                        r1 = lig3D.GetAtom(atom0).coords()
                        u = vecdiff(r1,mcoords)
                        dtheta = 5
                        dmin = 0
                        totiters = 0 
                        while totiters < 20:
                            lig3D = rotate_around_axis(lig3D,r1,u,dtheta)
                            d0 = lig3D.mindist(core3D)
                            if (d0 > dmin):
                                lig3Db = copy.deepcopy(lig3D)
                                dmin = d0
                            totiters += 1
                        lig3D = lig3Db
                    # get distance from bonds table or vdw radii
                    key = (metal,lig3D.GetAtom(atom0).sym+str(denticity))
                    if (key in MLbonds.keys()):
                        bondl = float(MLbonds[key][6-coord])
                    else:
                        bondl = m3D.GetAtom(0).rad + lig3D.GetAtom(atom0).rad
                    # get correct distance for center of mass
                    cmdist = bondl - distance(r1,mcoords)+distance(lig3D.centermass(),mcoords)
                    lig3D=setdistance(lig3D, mcoords, cmdist)
                elif (denticity == 2):
                    # connection atom
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(m3D,lig3D.GetAtom(atom0),m3D.GetAtom(totlig+1))
                    r1 = lig3D.GetAtom(atom0).coords()
                    # align center of mass to the middle
                    r21 = [a-b for a,b in zip(lig3D.GetAtom(catoms[1]).coords(),r1)]
                    r21n = [a-b for a,b in zip(m3D.GetAtom(totlig+2).coords(),r1)]
                    theta = 180*np.arccos(np.dot(r21,r21n)/(la.norm(r21)*la.norm(r21n)))/np.pi
                    u = np.cross(r21,r21n)
                    rrot = r1
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,rrot,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,rrot,u,theta-180)
                    d1 = distance(lig3D.GetAtom(catoms[1]).coords(),m3D.GetAtom(totlig+2).coords())
                    d2 = distance(lig3Db.GetAtom(catoms[1]).coords(),m3D.GetAtom(totlig+2).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one                 
                    # align center of mass
                    rm = [0.5*(a+b) for a,b in zip(lig3D.GetAtom(catoms[1]).coords(),r1)]
                    theta,u = rotation_params(r0,rm,lig3D.centermass())
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,rm,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,rm,u,theta-180)
                    d1 = lig3D.mindist(core3D)
                    d2 = lig3Db.mindist(core3D)
                    lig3D = lig3D if (d1 > d2)  else lig3Db # pick best one
                    r21 = vecdiff(r1,mcoords)
                    r21n = vecdiff(rm,mcoords)
                    costhb = np.dot(r21,r21n)/(la.norm(r21)*la.norm(r21n))+0.026
                    # get distance from bonds table or vdw radii
                    key = (metal,lig3D.GetAtom(atom0).sym+str(denticity))
                    if (key in MLbonds.keys()):
                        bondl = float(MLbonds[key][6-coord])
                    else:
                        bondl = m3D.GetAtom(0).rad + lig3D.GetAtom(atom0).rad
                    dbtotranslate = bondl*costhb + distance(rm,lig3D.centermass())
                    lig3D=setdistance(lig3D, mcoords, dbtotranslate)
                elif (denticity == 3):
                    atom0 = catoms[1]
                    baseatom = totlig + 2
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(m3D,lig3D.GetAtom(atom0),m3D.GetAtom(baseatom))
                    r1 = lig3D.GetAtom(atom0).coords()
                    r2 = lig3D.centermass() # center of mass
                    rrot = r1
                    theta,u = rotation_params(r0,r1,r2) 
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,rrot,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,rrot,u,theta-180)
                    d1 = distance(lig3D.centermass(),mcoords)
                    d2 = distance(lig3Db.GetAtom(atom0).coords(),mcoords)
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # rotate around secondary axis
                    u = vecdiff(lig3D.centermass(),mcoords)
                    r21 = vecdiff(lig3D.GetAtom(catoms[2]).coords(),mcoords)
                    alatom = baseatom+1 if totlig < 3 else baseatom-1 # alignment reference atom
                    r21n = vecdiff(m3D.GetAtom(alatom).coords(),mcoords)
                    theta = 180*np.arccos(np.dot(r21,r21n)/(la.norm(r21)*la.norm(r21n)))/np.pi
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,mcoords,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,mcoords,u,theta-180)
                    d1 = distance(lig3D.GetAtom(catoms[2]).coords(),m3D.GetAtom(baseatom-1).coords())
                    d2 = distance(lig3Db.GetAtom(catoms[2]).coords(),m3D.GetAtom(baseatom-1).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one                    
                    # get distance from bonds table or vdw radii
                    key = (metal,lig3D.GetAtom(atom0).sym+str(denticity))
                    if (key in MLbonds.keys()):
                        bondl = float(MLbonds[key][6-coord])
                    else:
                        bondl = m3D.GetAtom(0).rad + lig3D.GetAtom(atom0).rad
                    # get correct distance for center of mass
                    cmdist = bondl - distance(r1,mcoords)+distance(lig3D.centermass(),mcoords)
                    lig3D=setdistance(lig3D, mcoords, cmdist)
                # combine molecules
                core3D = core3D.combine(lig3D)
                del lig3D
            totlig += denticity
    return core3D

def structgen(installdir,args,rootdir,ligands,ligoc): 
    # get global variables class
    globs = globalvars()
    ############ LOAD DICTIONARIES ############
    mcores = readdict(installdir+'/Cores/cores.dict')
    licores = readdict(installdir+'/Ligands/ligands.dict')
    ancores = readdict(installdir+'/Anions/anions.dict')
    MLbonds = loaddata(installdir+'/Data/ML.dat')
    ########## END LOAD DICTIONARIES ##########
    strfiles = []
    # load molecule core
    core = core_load(installdir,args.core,mcores) 
    catoms = core.cat
    cgmolH=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=False)    
    cgmol=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=True)    
    initcore3D = convert_to_mol3D(core)
    if (ligands):
        # check if metal complex or not
        atno = cgmol.GetAtomWithIdx(catoms[0]).GetAtomicNum()
        if (args.core in globs.metals()): 
            core3D = mcomplex(core,ligands,ligoc,MLbonds,installdir,licores)
        else:
            # load characteristic group
            cgmolH=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=False)    
            cgmol=Chem.MolFromMolFile(installdir+'Cores/'+core.cg,removeHs=True)    
            totnumligands = 0 
            for i,ligand in enumerate(ligands):
                # load ligand with H2
                denticity = len(licores[ligand][1:])
                if denticity < 2:
                    lig = lig_load(installdir,ligand,licores)
                    # get occupancy
                    occ = ligoc[i] if i < len(ligoc) else 1
                    for j in range(0,int(occ)):
                        # get connection atom(s)
                        catom = catoms[totnumligands%len(catoms)]
                        # functionalize core
                        [core,cgmolH] = add_lig(core,catom,lig,cgmol,cgmolH)
                        totnumligands += 1 
            # convert molecule to mol3D
            core3D = convert_to_mol3D(core)
    else:
        # convert molecule to mol3D
        core3D = convert_to_mol3D(core)  
    # generate multiple geometric arrangements
    Nogeom = int(args.anionsnum) if args.anionsnum else 1 # number of different combinations
    ligname = ''.join("%s" % l[0:2] for l in ligands) # folder name
    if (args.anion):
        # load anion, add hydrogens and convert to mol3D
        anion = anion_load(installdir,args.anion,ancores)
        #anionHs = prepmol(anion)
        AllChem.EmbedMolecule(anion)
        an3D = convert_to_mol3D(anion)
        # get core size
        mindist = initcore3D.size
        # assign reference point
        Rp = core3D.centermass()
        # Generate base case (separated structures)
        base3D = protate(copy.deepcopy(an3D),Rp,[20*mindist,0.0,0.0])
        mols = []
        maxdist = mindist+float(args.maxd) # Angstrom, distance of non-interaction    
        mindist = mindist+float(args.mind) # Angstrom, distance of non-interaction
        fname = rootdir+'/'+args.core[0:3]+ligname+args.anion[0:2] 
        for i in range(0,Nogeom):        
            # generate random sequence of parameters for rotate()
            while True:
                R = 20+mindist if i==0 else random.uniform(mindist,maxdist) # get random distance, separated for i=0
                phi = random.uniform(0.0,360.0)
                theta = random.uniform(0.0,180.0)
                thetax = random.uniform(0.0,360.0)
                thetay = random.uniform(0.0,360.0)
                thetaz = random.uniform(0.0,360.0)
                # translate
                tr3D = protate(copy.deepcopy(an3D), Rp,[R,theta,phi])
                # rotate center of mass
                newmol = cmrotate(tr3D,[thetax,thetay,thetaz])
                # check for overlapping
                if not(newmol.overlapcheck(core3D,1)):
                    break
            # write new xyz file
            newmol.writemxyz(core3D,fname+str(i))
            # append filename
            strfiles.append(fname+str(i))
    else:
        fname = rootdir+'/'+args.core[0:3]+ligname
        core3D.writexyz(fname)
        strfiles.append(fname)
    print 'Generated ',Nogeom,' structures!'
    return strfiles



