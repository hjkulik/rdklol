#!/usr/bin/env python
''' 
Created on 05/13/15

@author: EI
'''
# Written by Tim Ioannidis May 13 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

####################################################
#########   Defines class of global    #############
########   variables that are shared   #############
##########    within the program       #############
####################################################

from math import sqrt 

# atoms dictionary contains atomic mass, atomic number, covalent radius, data from http://www.webelements.com/ (last accessed May 13th 2015)
amassdict = {'X':(1.0,0,0.77),'H':(1.0079,1,0.37),'C':(12.0107,6,0.77),'N':(14.0067,7,0.75),'O':(15.9994,8,0.73),'F':(18.9984,9,0.71),
        'P':(30.9738,15,1.06),'S':(32.065,16,1.02),'Cl':(35.453,17,0.99),'Ti':(47.867,22,1.36),'Cr':(51.9961,24,1.27),
        'Mn':(54.938,25,1.39),'Fe':(55.84526,26,1.25),'Ni':(58.4934,28,1.21),'Co':(58.9332,27,1.26),
        'Cu':(63.546,29,1.38),'Zn':(65.39,30,1.31),'Br':(79.904,35,1.14)}

metalslist = ['Sc','scandium','Ti','titanium','V','vanadium','Cr','chromium','Mn','manganese','Fe','iron','Co','cobalt','Ni','nickel',
            'Cu','copper','Zn','zinc','Y','yttrium','Zr','zirconium','Nb','niobium','Mo','molybdenum','Tc','technetium','Ru','ruthenium',
            'Rh','rhodium','Pd','palladium','Ag','silver','Cd','cadmium','Lu','lutetium','Hf','hafnium','Ta','tantalum',
            'W','tungsten','Re','rhenium','Os','osmium','Ir','iridium','Pt','platinum','Au','gold','Hg','mercury']

class globalvars:
    def amass(self):
        return amassdict
    def metals(self):
        return metalslist
