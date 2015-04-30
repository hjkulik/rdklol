#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script calculates the distance #############
########  between two molecules defined as   #############
########  the distance between their centers #############
########              of mass                #############
##########################################################

# import std modules
import argparse
import psycopg2
import urllib2
import glob
import os
import shutil
import time
import pickle
from math import sqrt
import sys
import copy
import random



# atoms dictionary contains atomic mass, atomic number, covalent radius
amass = {'H':(1.0079,1,0.37),'C':(12.0107,6,0.77),'N':(14.0067,7,0.75),'O':(15.9994,8,0.73),'F':(18.9984,9,0.71),
        'P':(30.9738,15,1.06),'S':(32.065,16,1.02),'Cl':(35.453,17,0.99),'Ti':(47.867,22,1.36),'Cr':(51.9961,24,1.27),
        'Mn':(54.938,25,1.39),'Fe':(55.84526,26,1.25),'Ni':(58.4934,28,1.21),'Co':(58.9332,27,1.26),
        'Cu':(63.546,29,1.38),'Zn':(65.39,30,1.31),'Br':(79.904,35,1.14)}

def distance(R1,R2):
    dx = R1[0] - R2[0] 
    dy = R1[1] - R2[1] 
    dz = R1[2] - R2[2] 
    d = sqrt(dx**2+dy**2+dz**2)
    return d

def centermass(mol):
    cm = [0.0,0.0,0.0]
    mass = 0.0
    for atom in mol:
        mass += amass[atom[0]][0]
        cm[0] += amass[atom[0]][0]*atom[1]
        cm[1] += amass[atom[0]][0]*atom[2]
        cm[2] += amass[atom[0]][0]*atom[3]
    cm[0] /= mass
    cm[1] /= mass
    cm[2] /= mass
    return cm

def calc_dist(xyzf,Nmol):    
    #load file and convert to molecule
    f = open(xyzf,'r')
    t=f.read()
    t = t.splitlines()
    mol0 = []
    mol1 = []
    NoAtoms = int(t[0])
    for i,line in enumerate(t[2:]):
        l = line.split(None)
        if (i<int(Nmol)):
            mol0.append([l[0],float(l[1]),float(l[2]),float(l[3])])
        else:
            mol1.append([l[0],float(l[1]),float(l[2]),float(l[3])])
    cm0 = centermass(mol0)
    cm1 = centermass(mol1)
    d=distance(cm0,cm1)
    return d
    
    
if __name__ == "__main__":
    if len(sys.argv) < 3 :
        print "\nUsage moldistance.py <file.xyz> <N> . Where N is the number of atoms for the 1st molecule\n"
        exit()
    xyzf = sys.argv[1]
    Nmol = sys.argv[2]
    d = calc_dist(xyzf,Nmol)
    print d
    
    
    
    
    
