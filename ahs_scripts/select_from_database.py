#!/bin/env python

import argparse
import pickle
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import SaltRemover
from rdkit.Chem import PandasTools
from rdkit.Chem import FunctionalGroups

import psycopg2
import pandas
import pandas.io.sql as psql

from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoDistMat
from rdkit.SimDivFilters import rdSimDivPickers


def query_database(args):
   
    limit=args.limit  

    #construct query string
    #query string is different for chembl_19 vs. emolecules
    if args.database=='chembl_19':
        squery = """SELECT molregno,m smiles,mol_numatoms(m),mol_send(m),bfp_to_binary_text(morganbv_fp(m)) 
                            FROM rdk.mols 
                            WHERE m @> {smarts}::qmol 
                            AND mol_numatoms(m) < {numatoms} 
                            LIMIT {lim}""".format(smarts='\''+smartsq+'\'',
                                                  numatoms=args.num_atoms,
                                                  lim=limit)
    elif args.database=='emolecules':
         squery = """SELECT id,m smiles,mol_numatoms(m),mol_send(m),bfp_to_binary_text(morganbv_fp(m)) 
                            FROM mols 
                            WHERE m @> {smarts}::qmol 
                            AND mol_numatoms(m) < {numatoms} 
                            LIMIT {lim}""".format(smarts='\''+smartsq+'\'',
                                                  numatoms=args.num_atoms,
                                                  lim=limit)
    else:
        print 'You have selected a database that this script does not understand.'

    print 'You are about to run the following query on {}:'.format(args.database)
    print '=======================================QUERY==========================================='
    print squery
    print '======================================================================================='
    #connect to database
    conn = psycopg2.connect(database=args.database,user='hjkulik',password='')

    #perform the database query into a Pandas dataframe
    conn.rollback()
    df = psql.read_sql(squery, con=conn)
    
    print 'Your query returned {} results'.format(df.shape[0])
    print '...'
    #fix the database objects
    df['mol_send']=df['mol_send'].map(lambda x: Chem.Mol(str(x)))
    df['mbv_fp']=df['bfp_to_binary_text'].map(lambda x: DataStructs.CreateFromBinaryText(str(x)))

    return df

def filter_ions(df):
    print 'You provided {} molecules to the salt remover'.format(df.shape[0])
    #Strip common ions out of molecule objects
    remover = SaltRemover.SaltRemover(defnData="[Li,Na,K,Rb,Cs,Mg,Ca,Sr,Ba,Zn,Cl,Br,F,I]")
    df['mol_strip'] = df['mol_send'].map(remover.StripMol)
    df['smilesf'] = df['mol_strip'].map(Chem.MolToSmiles)
    df['smilesf'] = df['smilesf'].map(lambda x: max(x.split('.'), key=len))
    print """CAUTION. You are removing ions and other fragments. However, the
            fingerprints used to calculate diversity were determined before removal.
            Consider recalculating fingerprints."""

    print 'Filter will try to remove duplicates after de-salting...'

    #remove duplicates remaining after removing counterions
    df.drop_duplicates(inplace=True,subset='smilesf')

    print 'After duplicate removal, there are {} molecules'.format(df.shape[0])
    print '...'
    return df

def filter_elements(df):
    print 'You provided {} molecules to the element filter'.format(df.shape[0])
    #create new database dropping all non HCNO elements
    #Could  do this in the first database call, 
    #but getting rid of chlorinated species and chlorine atoms makes this complicated
    if args.halogens:
        aquery = '[#3,#4,#5,#11,#15,#16,#32,#33,#34,#47,#51]'
        print 'You are allowing halogenated species to pass the filter.'
    else:
        aquery = '[#3,#4,#5,#9,#11,#15,#16,#17,#32,#33,#34,#35,#47,#51,#53]'
    
    dff = df.ix[~(df['mol_strip'] >= Chem.MolFromSmarts(aquery))]

    print 'After element filter, there are {} molecules'.format(dff.shape[0])
    print '...'
    return dff

def countSubstructMatches(mol, smartsq):
    matchlist = mol.GetSubstructMatches(Chem.MolFromSmarts(smartsq))
    return len(matchlist)

def single_match_filter(df):
    print 'You provided {} molecules to the single-match filter'.format(df.shape[0])
    #locate multiple matches to smarts pattern
    df['num_matches']=df['mol_strip'].map(lambda x: countSubstructMatches(x, smartsq))

    #create new frame containing only single matches
    dfm = df.ix[df['num_matches'] == 1]

    print 'You have {} molecules remaining with only a single match'.format(dfm.shape[0])
    print '...'
    return dfm

def remove_duplicates(df):
    print 'You provided {} molecules to the duplicate filter'.format(df.shape[0])
    #Now remove smarts pattern match, and create a smiles string for comparison
    #Will remove duplicate if carboxylic acid and carboxylate both in database
    patt = Chem.MolFromSmarts(smartsq)
    df['smiles_nopatt']=df['mol_strip'].map(lambda x: Chem.MolToSmiles(AllChem.DeleteSubstructs(x,patt)))

    df.drop_duplicates(inplace=True,subset='smiles_nopatt')

    print 'After removing duplicates, there are {} molecules remaining'.format(df.shape[0])
    print '...'
    return df

def pick_diverse_set(df,num_mols):
    print 'Resetting the index of the data'
    #reset the index of the data-frame
    dfr = df.reset_index()

    print 'Ther are {} molecules to calculate a distance matrix for'.format(dfr.shape[0])
    #Calculate Tanimoto Distance Matrix for the remaining molecules
    dm = GetTanimotoDistMat(dfr.mbv_fp.tolist())

    picker = rdSimDivPickers.MaxMinPicker()

    if num_mols >= dfr.shape[0]:
        print 'You are requesting more molecules than made it through the filters.'
        print 'Returning all the molecules.'

        return dfr
    else:
        #Should probably report some statistics here
        ids = picker.Pick(dm,dfr.shape[0],num_mols)

        dfo = dfr.ix[ids]

        print 'The diversity picker has selected {} molecules'.format(dfo.shape[0])
        print '...'
        return dfo
    

def output_results(df,fname):
    print 'Writing a pickled dictionary of molregnos and SMILES strings into file {}'.format(fname)
    #A crappy way to get some output
    outdict ={}
    if args.database=='chembl_19':
        for i in df.iterrows():
            outdict[i[1].molregno]=i[1].smilesf
    elif args.database=='emolecules':
        for i in df.iterrows():
            outdict[i[1].id]=i[1].smilesf
    else:
        print "I don't know how to output this database."
    
    with open(fname,'wb') as handle:
        pickle.dump(outdict, handle)

def output_grid(df,imagename):
    #use PandasTools.FrameToGridImage
    template = Chem.MolFromSmarts(smartsq)
    AllChem.Compute2DCoords(template)

    df['mol_strip'].apply(lambda x: AllChem.GenerateDepictionMatching2DStructure(x, template));
    df['matches']=df['mol_strip'].apply(lambda x: x.GetSubstructMatch(template))
    frameimage = PandasTools.FrameToGridImage(df, 
                                    column='mol_strip',
                                    legendsCol='molregno',
                                     molsPerRow=10,
                                     highlightAtomLists=list(df['matches']))
    frameimage.save(imagename)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--group", help="Functional group to search for in the database",default='CarboxylicAcid')
    parser.add_argument("-a", "--num_atoms", type=int,  
                    help="Maximum size of an molecule (in atoms) to be extracted from database",default=20)
    parser.add_argument("-c", "--halogens", action='store_true',  
                    help="Allow halogens to be included in selected ligands.",default=False)
    parser.add_argument("-n", "--num_mols", type=int, 
                    help="Number of divese molecules to extract from the full set",default=100)
    parser.add_argument("-d", "--database", 
                    help="Name of the local rdkit-enabled database to query (ex. chembl_19)",default='chembl_19')
    parser.add_argument("-f", "--outfile", 
                    help="File name for output list of molecules", default='output.pkl')
    parser.add_argument("-p", "--image", action='store_true', 
                    help="""Create a grid image of the molecules, up to 200 molecules.
                            Uses same file name as outfile, but with png extension""")
    parser.add_argument("-l", "--limit", type=int,
                    help="""Limits the number of molecules returned from the initial
                                            database query. Set to a low number to make the query fast 
                                            for debugging purposes.""", default=1000000) 
    args = parser.parse_args()
   
   #the functional groups are already coded into rdkit! TODO access these instead.
   # fgroup_smarts = {'carboxylate':'[CX3](=O)[O-]',
   #                  'acid':'[CX3](=O)[OX2H1]',
   #                  'both':'[CX3](=O)[OX1H0-,OX2H1]'}
   
   #Code below replaces orginal functional group queries:
    #CarboxylicAcid: C(=O)[O;H,-]
    #Should check if this is equivalent to old "both"
   
    fgs = FunctionalGroups.BuildFuncGroupHierarchy()
    from collections import namedtuple
    nt = namedtuple('pattern','smarts mol')
    def flattenFgs(fgs,res):
        if not fgs:
            return
        for x in fgs:
            res[x.label]=nt(x.smarts,x.pattern)
            flattenFgs(x.children,res)
    allFgDefs={}
    flattenFgs(fgs,allFgDefs)

    try:
        smartsq = allFgDefs[args.group].smarts
        print smartsq

    except KeyError,e :
        print """I can search for {}. 
                    You entered '{}'""".format(', '.join(allFgDefs.keys()), args.group)
        sys.exit(0) 

    df = query_database(args)
    df = filter_ions(df)
    df = filter_elements(df)
    df = single_match_filter(df)
    df = remove_duplicates(df)
    df = pick_diverse_set(df,args.num_mols)

    if args.image:
        imagename = args.outfile+'.png'
        output_grid(df.head(200),imagename)

    output_results(df,args.outfile)

