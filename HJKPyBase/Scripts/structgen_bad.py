''' 
Created on 12/20/14

@author: hjk
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalForceFields as FF
from rdkit.Chem import rdMolTransforms
from tcgen import mybash
import argparse
import psycopg2
import urllib2
import glob
import os
import shutil
import time
import pickle
import random
#This is the indium core geometry dictionary. It contains:
#-Name of PDB file containing In+carboxylate atoms.
#-Numbering of each carboxylate atom with C as the first atom, O/O as 2/3.
#-This is currently set to ensure ligands not overlap, could do conformer search for now.
incores={'bidentate':('bidentate_core.pdb',((1,3,6),(2,4,8),(0,5,7))),
         'monodentate':('smonodentate_core.pdb', ((1,3,6),(2,8,4),(0,7,5))),
         'two-bridge':('two_bridge_core.pdb', ((1,4,2),(0,3,5),(7,8,9),(11,13,16),(12,18,14),(10,15,17))),
         'three-bridge':('three_bridge_core.pdb', ((1,10,8),(2,9,6),(3,7,11),(0,5,4),(14,16,15),(17,18,19))), 
         'four-bridge':('four_bridge_core.pdb', ((1,7,5),(2,10,9),(14,15,16),(17,18,19),(3,8,11),(0,4,6)))}
basedir='/home/hjkulik/InStruct/Py/Cores/'
def core_check(userinput):
    """ Check if user input matches defined dictionary. Could add support for new user defined cores. Not for now.
    """
    if userinput not in incores:
       print("We didn't find the core structure: %s in the dictionary. Try again!\n" %(userinput))
       exit()
    if not glob.glob(basedir+incores[userinput][0]):
       print("We can't find the core structure file %s%s right now! Something is amiss.\n" %(basedir,incores[userinput][0]))
       exit()
    else:
       shutil.copy(basedir+incores[userinput][0],os.getcwd())
    return userinput

def lig_check(liginput,smiinput,dbinput):
    if smiinput:
       if 'C(=O)[O-]' not in smiinput[-9:]:
           print "Provide a carboxylate SMILES string that ends in C(=O)[O-]. Try again."
           exit()
       struct=smiinput
       if liginput:
          ligkey=liginput
       else:
          ligkey=smi_to_iupac(struct)
    elif liginput:
       ligkey=liginput
       if ligkey not in smilesdict:
           if not dbinput:
             print "Provide a name that matches a result in the SMILES dictionary or provide your own SMILES string. Try again."
             exit()
           else:
             struct='placeholder'
       else:
           struct=smilesdict[ligkey]
    return ligkey,struct

# Temporary smiles dictionary of handwritten SMILES strings.
# Want to put database support for chembl/emolecules in here shortly.
smilesdict={'benzoate': 'c1ccc(cc1)C(=O)[O-]',
            'myristate': 'CCCCCCCCCCCCCC(=O)[O-]',
            'acetate': 'CC(=O)[O-]'}

def smi_to_iupac(smi):
     """Tries to get structure name for a manual SMILES string."""
     try:
         url = 'http://cactus.nci.nih.gov/chemical/structure/'+smi+'/iupac_name';

         iupacName = urllib2.urlopen(url).read()
         return iupacName

     except urllib2.HTTPError, e:
         print "HTTP error: %d" % e.code
         return None
     except urllib2.URLError, e:
         print "Network error: %s" % e.reason.args[1]
         return None
     except:
         print "conversion failed for smiles "+ smi
         return None

def prepLig(lig):
    """Adds hydrogens to SMILES, turns it into a 3D, MMFF optimized object."""
    ligh=AllChem.AddHs(lig)
    AllChem.EmbedMolecule(ligh)
    AllChem.MMFFOptimizeMolecule(ligh)
    return ligh

def combAlignedOptLigHCore(core,lig,list):
    """Aligns ligand carboxylate to core carboxylate.
       Identifies which atom will need to be connected across ligand/core.
       Deletes ligand carboxylate.
       Combines ligand and core molecules to one molecule.
    """
    atomnums=lig.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]'))
    print "Alignment result: ", AllChem.AlignMol(lig,core,atomMap=zip(atomnums,list))
    connect_atom=lig.GetAtomWithIdx(lig.GetSubstructMatch(Chem.MolFromSmarts('*[CX3](=O)[OX1H0-,OX2H1]'))[0])
    connect_atom.SetProp('connect','Y')
    trunc=Chem.DeleteSubstructs(lig,Chem.MolFromSmarts('[CX3](=O)[OX2H1][H]'))
    if trunc.GetNumAtoms() == lig.GetNumAtoms():
       trunc=Chem.DeleteSubstructs(lig,Chem.MolFromSmarts('[CX3](=O)[OX1H0-]'))
    allatoms=trunc.GetAtoms()
    combo=Chem.CombineMols(core,trunc)
    return combo

def restoreBonds(combo,atom):
    """ Forms bonds between carboxylate neighbor atoms in ligand
        to carboxylate core. 
    """
    edcombo=AllChem.EditableMol(combo)
    for a in combo.GetAtoms():
        if a.HasProp('connect'):
           idx=a.GetIdx()
    edcombo.AddBond(atom,idx,Chem.rdchem.BondType.SINGLE) 
    newmol = edcombo.GetMol()
    return newmol

def deleteBonds(combo,atom):
    """ Remove bonds between carboxylate neighbor atoms and the metal ion.
        Not currently used.
    """
    edcombo=AllChem.EditablMol(combo)
    for a in combo.GetAtoms():
        if a.HasProp('disconnect'):
           idx=a.GetIdx()
    edcombo.RemoveBond(atom,idx)
    newmol = edcombo.GetMol()
    return newmol

def molQuality(mol):
    natom=mol.GetNumAtoms()
    conf=mol.GetConformer()
    minimum=1e100
    for i in range(0,natom-1):
        for j in range(i+1,natom):
         if mol.GetBondBetweenAtoms(i,j)==None:
            dist=rdMolTransforms.GetBondLength(conf,i,j)
            if dist<minimum:
                minimum=dist
    return minimum
if __name__ == "__main__":
    # Parse commandline arguments
    parser= argparse.ArgumentParser()
    parser.add_argument("-c","--core", help="core structure type")
    # want to 
    parser.add_argument("-l","--ligand", help="ligand structure name")
    parser.add_argument("-smi","--smiles", help="(optional) manual SMILES string")
    parser.add_argument("-q","--dbquery", help="(optional) database to query")
    parser.add_argument("-n","--dbnum", help="(optional) database list-based number to query", type=int)
    parser.add_argument("-rn","--regnum", help="(optional) database reg number to query", type=int)
    parser.add_argument("-f","--force", help="(optional) force build of structure with bad non-bonded distances", action="store_true")
    args=parser.parse_args()

    corekey=core_check(args.core)

    oldstructs=pickle.load(open('/home/hjkulik/InStruct/Py/Logs/setup.pkl','rb'))
    found=False
    numtorun=args.dbnum
    while not found: 
          if args.dbquery:
             # check for existence of pickle and read it in
             pickledir='/home/hjkulik/InStruct/Py/Pickles/'
             if glob.glob(pickledir+args.dbquery): 
                moldb=pickle.load(open(pickledir+args.dbquery,"rb"))
                moldblist=[]
                for key,value in moldb.iteritems():
                    moldblist.append([key,value]) 
                if numtorun >= 0:
                   ligkey,struct=str(moldblist[numtorun][0]),moldblist[numtorun][1]            
                   compound=corekey+'_'+ligkey
                   if compound not in oldstructs.keys(): 
                      found=True
                      oldstructs[compound]='Generated'
                   else:
                      found=False 
                      numtorun+=1
                      print("This job was already set up. Incrementing the db number up by 1 to %d!\n" %(numtorun))
                elif args.regnum:
                   ligkey=str(args.regnum)
                   struct=moldb[int(ligkey)]
                   compound=corekey+'_'+ligkey
                   if compound not in oldstructs.keys(): 
                      found=True
                      oldstructs[compound]='Generated'
                   else:
                      if not args.force:
                         found=False
                         print("This job was already set up. Exiting!\n")
                         exit()
                      else:
                         found=True
                         print struct
                         oldstructs[compound]='Generated'
                else: 
                   ligkey=str(random.choice(moldb.keys()))
                   struct=moldb[int(ligkey)]
                   compound=corekey+'_'+ligkey
                   if compound not in oldstructs.keys(): 
                      found=True
                      oldstructs[compound]='Generated'
                   else:
                      found=False
          else:
             ligkey,struct=lig_check(args.ligand,args.smiles,args.dbquery)
             compound=corekey+'_'+ligkey
             if compound not in oldstructs.keys(): 
                found=True
                oldstructs[compound]='Generated'
             else:
                found=True
                #print("This job was already set up. Exiting!\n")
                #exit()
 
    # Initializing single ligand model
       
    lig=Chem.MolFromSmiles(struct) 
    ligHs=prepLig(lig)
    liglist=[element for tupl in incores[corekey][1] for element in tupl]
    
    # Initializing core
    corekey=core_check(args.core)
    core=Chem.MolFromPDBFile(incores[corekey][0],sanitize=False)
    coreNat=len(core.GetAtoms())
    centerlist=[element for element in range(0,coreNat) if element not in liglist]
    coreatoms=core.GetAtoms()
    
    # Add and align ligands, restore bonds
    for matchset in range(0,len(incores[corekey][1])):
        list=incores[corekey][1][matchset]
        combo=combAlignedOptLigHCore(core,ligHs,list)
        core=restoreBonds(combo,list[0])
    
    # Generate PDB file
    pdbname=corekey+'_'+ligkey
    short=molQuality(core)
    if short < 1.5:
       if not args.force and (short < 1.0):
          print("Warning! Shortest non-bonded distance is %s. %s may be a bad structure! Not generated!\n" %(short,pdbname))
          exit()
       else:
          print("Warning! Shortest non-bonded distance is %s. %s may be a bad structure! Generated anyway!\n" %(short,pdbname))
    Chem.MolToPDBFile(core,pdbname+'.pdb',flavor=4)
    # Generate XYZ file ( has higher precision)
    mybash('babel -ipdb %s.pdb -oxyz %s.xyz --title "%s generated by structgen.py on %s"' %(pdbname,pdbname,pdbname,time.strftime("%m/%d/%y")))
    pickle.dump(oldstructs,open('/home/hjkulik/InStruct/Py/Logs/setup.pkl','wb'))
    ligtemp=open('ligstruct','w')
    ligtemp.write(ligkey) 
    ligtemp.close()
# Next tasks: querying from database, terachem script generation, terachem output collect

# Currently scrapped further FF-opt of complex. Would need more checks/balances for that.
# We need a main block and comments/guides to user...e.g the core/lig keys should be passed via some kind of input instead of hard coded.

