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
    
def mcomplex(core,ligs,ligoc,MLbonds,installdir,licores,ffopt):
    #################################################
    ####### functionalizes core with ligands ########
    ############## for metal complexes ##############
    #################################################
    coordbasef=['oct','tbp','thd','tri','li','one']
    metal = core.GetAtomWithIdx(0).GetSymbol()
    occs0 = []
    dentl = []
    toccs = 0 
    octa = False # flag for forced octahedral structures like porphyrines
    # find correct occurence for each ligand
    for i,lig in enumerate(ligs):
        dent_i = int(len(licores[lig][2:]))
        oc_i = int(ligoc[i]) if i < len(ligoc) else 1
        occs0.append(0)
        dentl.append(dent_i)
        for j in range(0,oc_i):
            if (toccs+dent_i <= 6):
                occs0[i] += 1
            toccs += dent_i
            if dent_i == 4:
                octa = True
    # sort by descending denticity (needed for adjacent connection atoms)
    indcs = [i[0] for i in sorted(enumerate(dentl), key=lambda x:x[1],reverse=True)]    
    ligands = [ligs[i] for i in indcs]
    occs = [occs0[i] for i in indcs]
    coord = min(toccs,6) # complex coordination
    coord = 6 if octa else coord
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
            denticity = len(licores[ligand][2:])
            if not(ligand=='x' or ligand =='X') and (totlig-1+denticity < coord):
                lig = lig_load(installdir,ligand,licores) # load ligand
                core3D.charge += lig.charge
                ###############################
                ### FORCE FIELD OPTIMIZATION ##
                if (ffopt and lig.mol):
                    AllChem.EmbedMolecule(lig)
                    if (ffopt == 'mmff'):
                        AllChem.MMFFOptimizeMolecule(lig)
                    else:
                        AllChem.UFFOptimizeMolecule(lig)
                ###############################
                if (lig.mol):
                    lig3D = convert_to_mol3D(lig) # get 3D ligand
                else:
                    lig3D = mol3D()
                    lig3D.readfromxyz(lig.f)
                catoms = lig.cat # connecting atom
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0 # initialize variables
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
                    lig3D=setcmdistance(lig3D, mcoords, cmdist)
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
                    lig3D=setcmdistance(lig3D, mcoords, dbtotranslate)
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
                    lig3D=setcmdistance(lig3D, mcoords, cmdist)
                elif (denticity == 4):
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(m3D,lig3D.GetAtom(atom0),m3D.GetAtom(1))
                    r0c = m3D.GetAtom(1).coords()
                    r1c = m3D.GetAtom(2).coords()
                    r2c = m3D.GetAtom(3).coords()
                    r0l = lig3D.GetAtom(catoms[0]).coords()
                    r1l = lig3D.GetAtom(catoms[1]).coords()
                    r2l = lig3D.GetAtom(catoms[2]).coords()
                    theta,uc = rotation_params(r0c,r1c,r2c) # normal vector to backbone plane
                    theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    lig3Db = copy.deepcopy(lig3D)
                    theta = 180*np.arccos(np.dot(uc,ul)/(la.norm(uc)*la.norm(ul)))/np.pi
                    u = np.cross(uc,ul)
                    # rotate around axis to match planes
                    lig3D = rotate_around_axis(lig3D,r0l,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,r0l,u,theta-180)
                    d1 = distance(lig3D.centermass(),mcoords)
                    d2 = distance(lig3Db.centermass(),mcoords)
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # rotate around secondary axis to match atoms
                    r0l = lig3D.GetAtom(catoms[0]).coords()
                    r1l = lig3D.GetAtom(catoms[1]).coords()
                    r2l = lig3D.GetAtom(catoms[2]).coords()
                    theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    rm = lig3D.centermass()
                    r1 = vecdiff(lig3D.centermass(),r0l)
                    r2 = vecdiff(mcoords,r0l)
                    theta = 180*np.arccos(np.dot(r1,r2)/(la.norm(r1)*la.norm(r2)))/np.pi
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,r0l,ul,theta)
                    lig3Db = rotate_around_axis(lig3Db,r0l,ul,theta-90)
                    d1 = distance(lig3D.centermass(),mcoords)
                    d2 = distance(lig3Db.centermass(),mcoords)      
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # translate to center of mass
                    dcm = vecdiff(mcoords,lig3D.centermass())
                    lig3D.translate(dcm)
                elif (denticity == 5):
                    # get center of mass 
                    ligc = mol3D()
                    for i in range(0,4): #5 is the non-planar atom
                        ligc.addatom(lig3D.GetAtom(catoms[i]))
                    # translate ligand to the middle of octahedral
                    lig3D.translate(vecdiff(mcoords,ligc.centermass()))
                    # get plane
                    r0c = m3D.GetAtom(totlig+1).coords()
                    r2c = m3D.GetAtom(totlig+2).coords()
                    r1c = mcoords
                    r0l = lig3D.GetAtom(catoms[0]).coords()
                    r2l = lig3D.GetAtom(catoms[1]).coords()
                    r1l = mcoords
                    theta,uc = rotation_params(r0c,r1c,r2c) # normal vector to backbone plane
                    theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    theta = 180*np.arccos(np.dot(uc,ul)/(la.norm(uc)*la.norm(ul)))/np.pi
                    u = np.cross(uc,ul)
                    lig3Db = copy.deepcopy(lig3D)
                    # rotate around axis to match planes
                    lig3D = rotate_around_axis(lig3D,mcoords,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,mcoords,u,180+theta)
                    d1 = distance(lig3D.GetAtom(catoms[4]).coords(),m3D.GetAtom(m3D.natoms-1+totlig).coords())
                    d2 = distance(lig3Db.GetAtom(catoms[4]).coords(),m3D.GetAtom(m3D.natoms-1+totlig).coords())
                    lig3D = lig3D if (d2 < d1)  else lig3Db # pick best one
                    # rotate around center axis to match backbone atoms
                    r0l = vecdiff(lig3D.GetAtom(catoms[0]).coords(),mcoords)
                    r1l = vecdiff(m3D.GetAtom(totlig+1).coords(),mcoords)
                    u = np.cross(r0l,r1l)
                    theta = 180*np.arccos(np.dot(r0l,r1l)/(la.norm(r0l)*la.norm(r1l)))/np.pi
                    lig3Db = copy.deepcopy(lig3D)
                    lig3D = rotate_around_axis(lig3D,mcoords,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,mcoords,u,theta-90)
                    d1 = distance(lig3D.GetAtom(catoms[0]).coords(),m3D.GetAtom(totlig+1).coords())
                    d2 = distance(lig3Db.GetAtom(catoms[0]).coords(),m3D.GetAtom(totlig+1).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                elif (denticity == 6):
                    # get center of mass 
                    ligc = mol3D()
                    for i in range(0,6):
                        ligc.addatom(lig3D.GetAtom(catoms[i]))
                    # translate metal to the middle of octahedral
                    core3D.translate(vecdiff(ligc.centermass(),mcoords))
                # cobine molecules
                core3D = core3D.combine(lig3D)
                del lig3D
            totlig += denticity
    return core3D

def ligadd(core,ligands,ligoc,installdir,licores,ffopt):
    #################################################
    ####### functionalizes core with ligands ########
    ################ for ferrocene ##################
    #################################################
    # load base core for coordination
    catoms = core.cat
    core3D = convert_to_mol3D(core)
    totlig = 0 
    Hlist = [] # list of hydrogens to be removed
    maxcoord = len(catoms) # maximum connected ligands
    mcoords = core3D.GetAtom(0).coords()
    nats = core3D.natoms
    for i,ligand in enumerate(ligands):
        lig = lig_load(installdir,ligand,licores) # load ligand
        denticity = lig.denticity
        if denticity < 2:
            # get occupancy
            occ = ligoc[i] if i < len(ligoc) else 1
            for j in range(0,int(occ)):
                if not(ligand=='x' or ligand =='X') and (totlig < maxcoord):
                    lig = lig_load(installdir,ligand,licores) # load ligand
                    core3D.charge += lig.charge
                    ###############################
                    ### FORCE FIELD OPTIMIZATION ##
                    if (ffopt):
                        AllChem.EmbedMolecule(lig)
                        if (ffopt == 'mmff'):
                            AllChem.MMFFOptimizeMolecule(lig)
                        else:
                            AllChem.UFFOptimizeMolecule(lig)
                    ###############################
                    if (lig.mol):
                        lig3D = convert_to_mol3D(lig) # get 3D ligand
                    else:
                        lig3D = mol3D()
                        lig3D.readfromxyz(lig.f)
                    # get connection atom(s)
                    catom = catoms[totlig]
                    lcatom = lig.cat[0] # connecting atom for ligand
                    lcatom3D = lig3D.GetAtom(lcatom)
                    catom3D = core3D.GetAtom(catom)
                    # get right distance
                    db = distance(mcoords,catom3D.coords())+catom3D.rad+lcatom3D.rad
                    # align to right axis
                    uc = vecdiff(catom3D.coords(),mcoords)
                    lig3D = aligntoaxis(lig3D,lcatom3D.coords(),mcoords,uc,db)
                    # align center of mass
                    r0 = lig3D.GetAtom(lcatom).coords()
                    u0 = vecdiff(lig3D.centermass(),r0)
                    theta = 180*np.arccos(np.dot(u0,uc)/(la.norm(u0)*la.norm(uc)))/np.pi
                    u = np.cross(u0,uc)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,r0,u,theta)
                    # cobine molecules
                    core3D = core3D.combine(lig3D)
                    # hydrogen to be removed
                    Hlist.append(core3D.getHsbyIndex(catom)[0])
                    totlig += 1
                else:
                    totlig += 1
        else:
            print "Multidentate ligands not available for core: "+args.core+". Skipping ligand "+ligand
    # remove extra hydrogens
    core3D.deleteatoms(Hlist)
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
    ########## START FUNCTIONALIZING ##########
    # load molecule core
    core = core_load(installdir,args.core,mcores)
    # get ox state/charge 
    initcore3D = convert_to_mol3D(core)
    if (ligands):
        # check if metal complex or not
        if (args.core in globs.metals()): 
            core3D = mcomplex(core,ligands,ligoc,MLbonds,installdir,licores,args.ff)
        else:
            # load characteristic group
            core3D = ligadd(core,ligands,ligoc,installdir,licores,args.ff)
    else:
        # convert molecule to mol3D
        core3D = convert_to_mol3D(core)
    # get charge 
    if (args.ox):
        core3D.charge += args.ox
    ############ END FUNCTIONALIZING ###########
    # generate multiple geometric arrangements
    Nogeom = int(args.anionsnum) if args.anionsnum and args.anion else 1 # number of different combinations
    ligname = ''.join("%s" % l[0:2] for l in ligands) # folder name
    if (args.anion):
        # load anion, add hydrogens and convert to mol3D
        anion = anion_load(installdir,args.anion,ancores)
        #anionHs = prepmol(anion)
        AllChem.EmbedMolecule(anion)
        AllChem.MMFFOptimizeMolecule(anion)
        an3D = convert_to_mol3D(anion)
        # get core size
        mindist = initcore3D.size
        #mindist = core3D.size
        # assign reference point
        Rp = core3D.centermass()
        # Generate base case (separated structures)
        base3D = protate(copy.deepcopy(an3D),Rp,[20*mindist,0.0,0.0])
        mols = []
        maxdist = mindist+float(args.maxd) # Angstrom, distance of non-interaction    
        mindist = mindist+float(args.mind) # Angstrom, distance of non-interaction
        fname = rootdir+'/'+args.core[0:3]+ligname+args.anion[0:2] 
        core3D.charge += anion.charge
        for i in range(0,Nogeom):        
            # generate random sequence of parameters for rotate()
            totits = 0 
            while True:
                R = random.uniform(mindist,maxdist) # get random distance, separated for i=0
                phi    = random.uniform(0.0,360.0)
                rann = random.uniform(0.0,180.0)
                theta  =  rann*np.exp((rann-180.0)/100.0) # biased towards 0 
                #theta = random.uniform(0.0,180.0)
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
                if totits > 200:
                    print "WARNING: Overlapping in molecules for file "+fname+str(i)
                    break 
                totits += 1
            # write new xyz file
            newmol.writemxyz(core3D,fname+str(i))
            # append filename
            strfiles.append(fname+str(i))
    else:
        fname = rootdir+'/'+args.core[0:3]+ligname
        core3D.writexyz(fname)
        strfiles.append(fname)
    print 'Generated ',Nogeom,' structures!'
    return strfiles,core3D.charge



