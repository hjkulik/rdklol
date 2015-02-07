''' 
Created on 12/27/14

@author: hjk
'''
import argparse
import glob
import sys
import os
import subprocess
from tcgen import mybash,tc_check,job_check
from structgen import lig_check,core_check

if __name__ == "__main__":
    # Parse commandline arguments
    parser= argparse.ArgumentParser()
    #COME BACK: this needs to be replaced with stuff that pertains to structgen - a ligand and a core/which ligand/core sets.
    # imported over from structgen:
    parser.add_argument("-c","--core", help="core structure type")
    parser.add_argument("-l","--ligand", help="ligand structure name")
    parser.add_argument("-smi","--smiles", help="(optional) manual SMILES string")
    parser.add_argument("-d","--dbquery", help="(optional) database to query")
    # need to make these seamless - ie jobname is core_ligand and ligand needs to be named from smiles/dbquery
    parser.add_argument("-mem","--memory", help="memory reserved per thread for job file (default: 2G.")
    parser.add_argument("-t","--time", help="time requested for queueing system (default: 168hrs).")
    parser.add_argument("-q","--queue", help="queue wildcards (default: [dg]*).")
    parser.add_argument("-g","--gpus", help="number of GPUS (default: 1).")
    parser.add_argument("-m","--method", help="electronic structure approach for terachem job. Specify UHF/UDFT to enforce levelshifting and unrestricted calculation for OS singlets (default: b3lyp).")
    parser.add_argument("-b","--basis", help="basis for terachem job (default: LACVP*).")
    parser.add_argument("-s","--spin", help="spin for system (default: singlet).")
    parser.add_argument("-ch","--charge", help="charge for system (default: neutral).")
    parser.add_argument("-x","--extra", nargs="+", help="extra arguments in syntax tckeyword=value. suggestions: nstep (opt), maxit (scf), scf=diis+a, min_tolerance>4.5e-4/min_tolerance_e>1.0e-6, min_coordinates=cartesian, orbitalswrtfrq=<num>,multibasis=<basis file>.")
    parser.add_argument("-tc","--terachem", help="full path to custom terachem installation.")
    # TBD, could do groupings of smart keywords at some point - eg 'floppy' or 'mustconverge' but not for now.
    args=parser.parse_args()
    # Initialize the jobscrparams dictionary with mandatory/useful keywords.
    jobscrparams={'memory': '2G',
                  'time': '168:00:00',
                  'queue': '[dg]*',
                  'gpus': '1'}
    # Set queue wildcard based on CLI.
    if args.queue:
       qresult=mybash('qstat -f -q %s\n' %(args.queue))
       if not qresult:
          qresult_wc=mybash('qstat -f -q %s*\n' %(args.queue))
          if not qresult_wc: 
             print("Your queue pattern/wildcard: %s was invalid. Try again." %(args.queue))
             exit()
          else: 
             jobscrparams['queue']=args.queue+'*'
       else:
             jobscrparams['queue']=args.queue
    if args.extra:
       for elem in range(0,len(args.extra)):
          key,val=args.extra[elem].split('=')
          jobscrparams[key]=val
    corekey=core_check(args.core)
    ligkey,struct=lig_check(args.ligand,args.smiles,args.dbquery) 
    jobname=corekey+'_'+ligkey
    jobname,jobtrunc,jobfile=job_check(jobname)
    ld_extra,tc_exe,modname=tc_check(args.terachem)
    output=open('qscr','w')
    output.write('#$ -S /bin/bash\n')
    output.write('#$ -N %s\n' %(jobtrunc))
    output.write('#$ -l gpus=1\n')
    output.write('#$ -cwd\n')
    output.write('#$ -R y\n')
    output.write('#$ -l h_rt=%s\n' %(jobscrparams['time']))
    output.write('#$ -l h_rss=%s\n' %(jobscrparams['memory']))      
    output.write('#$ -q %s\n' %(jobscrparams['queue']))
    output.write('#$ -pe smp %s\n' %(jobscrparams['gpus']))
    output.write('# -fin %s.*\n' %(jobname))
    if 'scrdir' in jobscrparams:
       scrdir=jobscrparams['scrdir']
    else:
       scrdir='scr/'
    output.write('# -fout %s\n' %(scrdir))
    output.write('module load intel\n')
    #TBD: handle for specifying which cuda version
    output.write('module load cuda\n')
    if modname:
       output.write('module load %s\n' %(modname)) 
    if ld_extra:
       output.write('export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n' %(ld_extra))
    output.write('export OMP_NUM_THREADS=%s\n' %(jobscrparams['gpus']))
    output.write('%s %s.in > $SGE_O_WORKDIR/%s.out\n' %(tc_exe,jobname,jobname))
    #COME BACK: Include output parser/setup logger 
    #COME BACK: Need to tell it when to re-run the job, e.g if it fails and which params to change
    output.close()
