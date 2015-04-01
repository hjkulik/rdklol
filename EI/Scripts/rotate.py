#!/usr/bin/env python
import sys
import numpy as np
from numpy import linalg as LA
from math import pi ,sin, cos
from Classes.mol3D import mol3D
from Classes.atom3D import atom3D

# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script gets a set of atoms and #############
########   rotates/translates the molecule   #############
########  with respect to a point in space   #############
########     in order increase the radial    #############
########    distance between the center of   #############
######## mass of the molecule and the point  #############
######## and if requested rotates by polar   #############
######## angle theta and azimuthal angle phi #############
######## Furthermore, it rotates a molecule  #############
######## around its own center of mass if    #############
########   requested (x,y,z axes rotation)   #############
##########################################################

# INPUT for protate(Mol, Rr, Delta) or cmrotate(Mol, Delta)
# Mol: contains atom types and coordinates
# e.g.: mol =[['Fe',0.0,0.0,0.0],['C', -1.25700,       -0.18100,        1.59400],
#    ['C',         -0.57900,        1.06300,        1.63700]]
# Rr: coordinates of reference point for rotation/translation
# Delta: if len(args)==3 contains D->R, Dtheta, Dphi for rotation translation around reference point (distance R, angle theta, angle phi)
# Delta: if (len(args)==2 contains D->Dthetax, Dthetay, Dthetaz for rotation around center of mass
# e.g. protate(Mol,[0.0 0.0 0.0],[2 45.0 90.0] translates 2A around  [0.0 0.0 0.0] and rotates by theta=45.0 and phi=90.0 
# e.g cmrotate(Mol,[60.0 30.0 0.0]) rotates around center of mass and x-axis by thetax=60.0, y-axis thetay=30.0 and z thetaz=0.0

def translate_around_Rp(Mol,Pr,D):
    # get center of mass
    pmc = Mol.centermass()
    # get translation vector that corresponds to new coords
    Rt = PointTranslateSph(Pr,pmc,D)
    # translate molecule
    Mol.translate(Rt)
    return Mol
    
def rotate_around_cm(Mol,D):
    pmc = Mol.centermass()
    for atom in Mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(pmc,atom.coords(),D)
        atom.set_coords(Rt)
    return Mol

def PointTranslateSph(Pr,p0,D):    
    from math import cos, sin, sqrt 
    # translate to origin
    ps=[0,0,0]
    ps[0] = p0[0] - Pr[0]
    ps[1] = p0[1] - Pr[1]
    ps[2] = p0[2] - Pr[2]  
    # get current spherical coords
    r0 = LA.norm(ps) 
    if (r0 < 1e-16):
        theta0 = 0.5*pi
        phi0 = 0 
    else :
        theta0 = np.arccos(ps[2]/r0)
        phi0 = np.arctan2(ps[1],ps[0]) # atan doesn't work due to signs
    # get translation vector
    p = [0,0,0]
    p[0] = -D[0]*sin(theta0+D[1])*cos(phi0+D[2]) - p0[0]
    p[1] = -D[0]*sin(theta0+D[1])*sin(phi0+D[2]) - p0[1]
    p[2] = -D[0]*cos(theta0+D[1]) - p0[2]
    return p
    
def PointRotateSph(Pr,p0,D):    
    from math import cos, sin, sqrt 
    # translate to origin
    ps=[0,0,0]
    ps[0] = p0[0] - Pr[0]
    ps[1] = p0[1] - Pr[1]
    ps[2] = p0[2] - Pr[2]  
    # build 3D rotation matrices around x,y,z axes
    Mx=[[1, 0, 0],[0, cos(D[0]), -sin(D[0])],[0, sin(D[0]), cos(D[0])]]
    My=[[cos(D[1]), 0, sin(D[1])],[0, 1, 0],[-sin(D[1]), 0, cos(D[1])]]
    Mz=[[cos(D[2]), -sin(D[2]), 0],[sin(D[2]), cos(D[2]), 0],[0, 0, 1]]
    # get full rotation matrix
    M = np.array(np.mat(Mx)*np.mat(My)*np.mat(Mz))
    p=[0.0, 0.0, 0.0]
    # rotate atom and translate it back from origin
    p[0] = M[0][0]*ps[0] + M[0][1]*ps[1] + M[0][2]*ps[2] + Pr[0]
    p[1] = M[1][0]*ps[0] + M[1][1]*ps[1] + M[1][2]*ps[2] + Pr[1]
    p[2] = M[2][0]*ps[0] + M[2][1]*ps[1] + M[2][2]*ps[2] + Pr[2]
    return p

def protate(Mol, Rr, Delta):
    # rotate/translate around reference point
    # convert to rad
    Delta[0] = float(Delta[0])
    Delta[1] = (float(Delta[1])/180)*pi
    Delta[2] = (float(Delta[2])/180)*pi
    # rotate/translate around reference point
    Mol = translate_around_Rp(Mol,Rr,Delta)
    return Mol

def cmrotate(Mol, Delta):
    # rotate around center of mass    
    # convert to rad
    Delta[0] = (float(Delta[0])/180)*pi
    Delta[1] = (float(Delta[1])/180)*pi
    Delta[2] = (float(Delta[2])/180)*pi
    Mol = rotate_around_cm(Mol,Delta)
    return Mol

if __name__ == "__main__":
    if len(sys.argv) < 5 :
        print '\n Usage for translate/rotate around RP: python rotate.py input.xyz <Reference point> <Dr> <Dtheta> <Dphi> '
        print 'Example : python rotate.py input.xyz 0.0 0.0 0.0 0.5 30.0 60.0 \n'
        print 'OR\n'
        print '\n Usage for rotation around CM: python rotate.py input.xyz <Dthetax> <Dthetay> <Dthetaz>'
        print 'Example : python rotate.py input.xyz 20.0 30.0 40.0\n'
        sys.exit()
    filename = sys.argv[1]   
    # initialize mol
    mol = []
    with open(filename) as f:
        data = f.readlines()
    f.closed
    h = filename.split('.xyz')
    output=h[0]+'_rotated.xyz'
    f = open(output,'w')    
    f.write(str(len(data)-2)+'\n\n')   
    for idx,line in enumerate(data):
        newline = " ".join(line.split())
        fin = newline.split(' ')
        if idx > 1 :
           # reference point
           R = [float(fin[1]),float(fin[2]),float(fin[3])]        
           mol.append([fin[0],R[0],R[1],R[2]])
    if len(sys.argv)==8:
        Rr=[float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])]
        D=[float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])]
        Mol = protate(mol, Rr, D)    
    elif len(sys.argv)==5:
        D=[float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])]
        Mol = cmrotate(mol, D)    
    else:
        print 'Wrong number of input arguments. Exiting..'
        sys.exit()
    for atom in Mol:
        f.write(atom[0] + '\t' + str(atom[1]) + '\t'+ str(atom[2]) + '\t'+str(atom[3])+'\t' + '\n')  
    f.closed    
