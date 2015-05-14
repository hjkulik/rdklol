#!/usr/bin/env python
''' 
Created on 12/20/14

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
########   Defines class of 3D atoms that    #############
########     will be used to manipulate      #############
########      coordinates within RDkit       #############
##########################################################

from math import sqrt 
from Classes.globalvars import globalvars


class atom3D:
    """ Class atom3D represents an atom with its coordinates for
    easy manipulation in 3D space """
    def __init__(self,Sym,xyz):
        """ Create a new atom object """
        self.sym = Sym
        globs = globalvars()
        amass = globs.amass()
        if Sym not in amass:
            print("We didn't find the atomic mass of %s in the dictionary. Assigning default value of 12!\n" %(Sym))
            self.mass = 12 # default atomic mass
            self.atno = 6 # default atomic number
            self.rad = 0.75 # default atomic radius
        else:
            self.mass = amass[Sym][0] # atomic mass
            self.atno = amass[Sym][1] # atomic number
            self.rad = amass[Sym][2] # atomic covalent radius
        self.__xyz = xyz # coords
    def symbol(self):   # a class method
        return self.sym
    def ismetal(self):
        if self.sym in metals:
            return True
        else:
            return False
    def coords(self):   # a class method
        x,y,z = self.__xyz
        return [x,y,z]
    def set_coords(self,xyz):
        self.__xyz[0] = xyz[0]
        self.__xyz[1] = xyz[1]
        self.__xyz[2] = xyz[2]
    def distance(self,atom2):
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return sqrt(dx*dx+dy*dy+dz*dz)
    def distancev(self,atom2):
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return [dx,dy,dz]
    def translate(self,dxyz):   # a class method
        x,y,z = self.__xyz
        self.__xyz[0] = x + dxyz[0]
        self.__xyz[1] = y + dxyz[1]
        self.__xyz[2] = z + dxyz[2]
    def __repr__(self):
        """ when calls mol3D object without attribute e.g. t """
        ss = "\nClass atom3D has the following methods:\n"
        for method in dir(self):
            if callable(getattr(self, method)):
                ss += method +'\n'
        return ss
