from math import sqrt
import sys
at1=sys.argv[2]
at2=sys.argv[3]
myfile=open(sys.argv[1],'r').readlines()
numAtoms=int(myfile[0].strip('\n'))
at1pos=len(myfile)-numAtoms+int(at1)-1
at2pos=len(myfile)-numAtoms+int(at2)-1
x1,y1,z1=map(float,myfile[at1pos].strip('\n').split()[1:4])
x2,y2,z2=map(float,myfile[at2pos].strip('\n').split()[1:4])
dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
print dist
