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

# atoms dictionary contains atomic mass, atomic number, covalent radius (URL!!!!)
amass = {'X':(1.0,0,0.77),'H':(1.0079,1,0.37),'C':(12.0107,6,0.77),'N':(14.0067,7,0.75),'O':(15.9994,8,0.73),'F':(18.9984,9,0.71),
        'P':(30.9738,15,1.06),'S':(32.065,16,1.02),'Cl':(35.453,17,0.99),'Ti':(47.867,22,1.36),'Cr':(51.9961,24,1.27),
        'Mn':(54.938,25,1.39),'Fe':(55.84526,26,1.25),'Ni':(58.4934,28,1.21),'Co':(58.9332,27,1.26),
        'Cu':(63.546,29,1.38),'Zn':(65.39,30,1.31),'Br':(79.904,35,1.14)}

class atom3D:
    """ Class atom3D represents an atom with its coordinates for
    easy manipulation in 3D space """
    def __init__(self,Sym,xyz):
        """ Create a new atom object """
        self.sym = Sym
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
