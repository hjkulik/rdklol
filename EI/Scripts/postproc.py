#!/usr/bin/env python
''' 
Created on 07/07/2015

@author: EI
'''
# Written by Tim Ioannidis Mar 12 2015 for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
########  Medium level script that coordinates  ##########
############    post-processing of results   #############
##########################################################

from Scripts.structgen import *
from Scripts.tcgen import tcgen
from Scripts.jobgen import jobgen
from Scripts.io import *
import argparse
import sys
import os
import shutil
from collections import Counter
import itertools

def postproc()

