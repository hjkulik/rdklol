#!/usr/bin/env python
import sys
import copy
import numpy as np
from numpy import linalg as LA
from math import pi ,sin, cos, sqrt
from Classes.mol3D import mol3D
from Classes.atom3D import atom3D

# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
############# This script performs a lot of ##############
######### different operations and manipulates ###########
################### atoms in 3D space ####################
##########################################################

def ReflectPlane(u,r,Rp):
    #################################################
    ####### contains roto-translation matrix ########
    ####### for reflection of r through plane #######
    ########## with normal u and point rp ###########
    #################################################
    # INPUT
    #   - u: normal vector to plane [ux,uy,uz]
    #   - Rp: reference point [x,y,z]
    #   - r: point to be reflected [x,y,z]
    # OUTPUT
    #   - rn: reflected point [x,y,z]
    un = LA.norm(u)
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    # construct augmented vector rr = [r;1]
    d = -u[0]*Rp[0]-u[1]*Rp[1]-u[2]*Rp[2]
    # reflection matrix
    R=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    rn = [0,0,0]
    R[0][0] = 1-2*u[0]*u[0]
    R[0][1] = -2*u[0]*u[1] 
    R[0][2] = -2*u[0]*u[2] 
    R[0][3] = -2*u[0]*d
    R[1][0] = -2*u[1]*u[0] 
    R[1][1] = 1-2*u[1]*u[1] 
    R[1][2] = -2*u[1]*u[2] 
    R[1][3] = -2*u[1]*d
    R[2][0] = -2*u[2]*u[0]
    R[2][1] = -2*u[1]*u[2]
    R[2][2] = 1-2*u[2]*u[2] 
    R[2][3] = -2*u[2]*d
    R[3][3] = 1 
    # get new point
    rn[0] = R[0][0]*r[0]+R[0][1]*r[1]+R[0][2]*r[2] + R[0][3] 
    rn[1] = R[1][0]*r[0]+R[1][1]*r[1]+R[1][2]*r[2] + R[1][3] 
    rn[2] = R[2][0]*r[0]+R[2][1]*r[1]+R[2][2]*r[2] + R[2][3] 
    return rn

def PointRotateAxis(u,rp,r,theta):
    #################################################
    ####### contains roto-translation matrix ########
    ######### for rotation of r about axis u ########
    ######## through point rp by angle theta ########
    #################################################
    # INPUT
    #   - u: direction vector of axis [ux,uy,uz]
    #   - rp: reference point [x,y,z]
    #   - r: point to be reflected [x,y,z]
    #   - theta: angle of rotation in degrees
    # OUTPUT
    #   - rn: reflected point [x,y,z]
    # construct augmented vector rr = [r;1]
    rr = r
    rr.append(1)
    # rotation matrix about arbitrary line through rp
    R=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    rn = [0,0,0]
    R[0][0] = cos(theta)+u[0]**2*(1-cos(theta))
    R[0][1] = u[0]*u[1]*(1-cos(theta))-u[2]*sin(theta)
    R[0][2] = u[0]*u[2]*(1-cos(theta))+u[1]*sin(theta)
    R[0][3] = (rp[0]*(u[1]**2+u[2]**2)-u[0]*(rp[1]*u[1]+rp[2]*u[2]))*(1-cos(theta))
    R[0][3] += (rp[1]*u[2]-rp[2]*u[1])*sin(theta)
    R[1][0] = u[1]*u[0]*(1-cos(theta))+u[2]*sin(theta)
    R[1][1] = cos(theta)+u[1]**2*(1-cos(theta))
    R[1][2] = u[1]*u[2]*(1-cos(theta))-u[0]*sin(theta)
    R[1][3] = (rp[1]*(u[0]**2+u[2]**2)-u[1]*(rp[0]*u[0]+rp[2]*u[2]))*(1-cos(theta))
    R[1][3] += (rp[2]*u[0]-rp[0]*u[2])*sin(theta)
    R[2][0] = u[2]*u[0]*(1-cos(theta))-u[1]*sin(theta)
    R[2][1] = u[2]*u[1]*(1-cos(theta))+u[0]*sin(theta)
    R[2][2] = cos(theta)+u[2]**2*(1-cos(theta))
    R[2][3] = (rp[2]*(u[0]**2+u[1]**2)-u[2]*(rp[0]*u[0]+rp[1]*u[1]))*(1-cos(theta))
    R[2][3] += (rp[0]*u[1]-rp[1]*u[0])*sin(theta)
    R[3][3] = 1
    # get new point
    rn[0] = R[0][0]*r[0]+R[0][1]*r[1]+R[0][2]*r[2] + R[0][3]
    rn[1] = R[1][0]*r[0]+R[1][1]*r[1]+R[1][2]*r[2] + R[1][3]
    rn[2] = R[2][0]*r[0]+R[2][1]*r[1]+R[2][2]*r[2] + R[2][3]
    return rn
    
def PointTranslateSph(Rp,p0,D):    
    ##################################################
    ####### converts to spherical coordinates ########
    ######### translates/rotates and converts ########
    ################ back to cartesian ###############
    ##################################################
    # INPUT
    #   - Rp: reference point [x,y,z]
    #   - p0: point to be translated [x,y,z]
    #   - D: [radial distance, polar theta, azimuthal phi]
    # OUTPUT
    #   - p: translated point [x,y,z] 
    # translate to origin
    ps=[0,0,0]
    ps[0] = p0[0] - Rp[0]
    ps[1] = p0[1] - Rp[1]
    ps[2] = p0[2] - Rp[2]  
    # get current spherical coords
    r0 = LA.norm(ps) 
    if (r0 < 1e-16):
        theta0 = 0.5*pi
        phi0 = 0 
    else :
        theta0 = np.arccos(ps[2]/r0) #changed sign TIM IOANNIDIS
        phi0 = np.arctan2(ps[1],ps[0]) # atan doesn't work due to signs
    # get translation vector
    p = [0,0,0]
    p[0] = -D[0]*np.sin(theta0+D[1])*np.cos(phi0+D[2]) - p0[0]
    p[1] = -D[0]*np.sin(theta0+D[1])*np.sin(phi0+D[2]) - p0[1]
    p[2] = -D[0]*np.cos(theta0+D[1]) - p0[2]
    return p
    
def PointRotateSph(Rp,p0,D):    
    #############################################
    ############# contains rotation #############
    ########## matrix about x,y,z axes ##########
    #############################################
    # INPUT
    #   - Rp: reference point [x,y,z]
    #   - p0: point to be rotated [x,y,z]
    #   - D: [theta-x, theta-y, theta-z]
    # OUTPUT
    #   - p: rotated point
    # translate to origin (reference)
    ps=[0,0,0]
    ps[0] = p0[0] - Rp[0]
    ps[1] = p0[1] - Rp[1]
    ps[2] = p0[2] - Rp[2]  
    # build 3D rotation matrices about x,y,z axes
    Mx=[[1, 0, 0],[0, cos(D[0]), -sin(D[0])],[0, sin(D[0]), cos(D[0])]]
    My=[[cos(D[1]), 0, sin(D[1])],[0, 1, 0],[-sin(D[1]), 0, cos(D[1])]]
    Mz=[[cos(D[2]), -sin(D[2]), 0],[sin(D[2]), cos(D[2]), 0],[0, 0, 1]]
    # get full rotation matrix
    M = np.array(np.mat(Mx)*np.mat(My)*np.mat(Mz))
    p=[0.0, 0.0, 0.0]
    # rotate atom and translate it back from origin
    p[0] = M[0][0]*ps[0] + M[0][1]*ps[1] + M[0][2]*ps[2] + Rp[0]
    p[1] = M[1][0]*ps[0] + M[1][1]*ps[1] + M[1][2]*ps[2] + Rp[1]
    p[2] = M[2][0]*ps[0] + M[2][1]*ps[1] + M[2][2]*ps[2] + Rp[2]
    return p
    
def reflect_through_plane(mol,u,Rp):
    ############################################
    ##### reflects molecule through plane ######
    ##### with normal u and a point rp #########
    ############################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - u: normal vector to plane
    #   - rp: reference point on plane
    # OUTPUT
    #   - mol: reflected molecule
    un = LA.norm(u)
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = ReflectPlane(u,atom.coords(),Rp)
        atom.set_coords(Rt)
    return mol

def rotate_around_axis(mol,Rp,u,theta):
    ###############################################
    ########## rotates molecule about #############
    ####### arbitrary axis with direction #########
    ############# u through point rp ##############
    ################ by angle theta ###############
    ###############################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - Rp: reference point [x,y,z]
    #   - u: direction vector [vx,vy,vz]
    #   - theta: angle to rotate in degrees
    # OUTPUT
    #   - mol: translated molecule
    un = LA.norm(u)
    theta = (theta/180.0)*pi
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateAxis(u,Rp,atom.coords(),theta)
        atom.set_coords(Rt)
    return mol

def setdistance(mol, Rp, bond):
    #################################################
    ########## sets distance of molecule ############
    ############## from reference point #############
    #################################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - Rp: reference point [x,y,z]
    #   - bond: final bond length between molecule, center of mass
    # OUTPUT
    #   - mol: translated molecule
    # get float bond length
    bl = float(bond)
    # get center of mass
    cm = mol.centermass()
    # get unit vector through line r = r0 + t*u
    u = [a-b for a,b in zip(cm,Rp)]
    t = bl/LA.norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for centermass
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-cm[0]
    dxyz[1] = Rp[1]+t*u[1]-cm[1]
    dxyz[2] = Rp[2]+t*u[2]-cm[2]
    for atom in mol.atoms:
        # Translate atoms
        atom.translate(dxyz)
    return mol

def protate(mol, Rr, D):
    ##############################################
    ######## translates/rotates molecule #########
    ########### about reference point ############
    ##############################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - Rr: reference point [x,y,z]
    #   - D: [radial distance, polar theta, azimuthal phi]
    # OUTPUT
    #   - mol: translated molecule
    # rotate/translate about reference point
    # convert to rad
    D[0] = float(D[0])
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # rotate/translate about reference point
    # get center of mass
    pmc = mol.centermass()
    # get translation vector that corresponds to new coords
    Rt = PointTranslateSph(Rr,pmc,D)
    # translate molecule
    mol.translate(Rt)
    return mol

def cmrotate(mol, D):
    ########################################
    ########## rotates molecule ############
    ######## about center of mass ##########
    ########################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - D: [theta-x, theta-y, theta-z]
    # OUTPUT
    #   - mol: translated molecule
    # convert to rad
    D[0] = (float(D[0])/180.0)*pi
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # perform rotation
    pmc = mol.centermass()
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(pmc,atom.coords(),D)
        atom.set_coords(Rt)
    return mol
    
def pmrotate(mol, Rp, D):
    ########################################
    ########## rotates molecule ############
    ######## about reference point #########
    ########################################
    # INPUT
    #   - mol: molecule to be manipulated
    #   - Rp: reference point [x,y,z]
    #   - D: [theta-x, theta-y, theta-z]
    # OUTPUT
    #   - mol: translated molecule
    # convert to rad
    D[0] = (float(D[0])/180.0)*pi
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # perform rotation
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(Rp,atom.coords(),D)
        atom.set_coords(Rt)
    return mol   
