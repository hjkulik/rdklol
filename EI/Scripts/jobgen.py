#!/usr/bin/env python
''' 
Created on 12/27/14

@author: hjk
'''
import argparse
import glob
import sys
import os
import subprocess
from tcgen import mybash
from tcgen import job_check
from tcgen import mybash
from tcgen import tc_check

def jobgen(args,jobsdir,lig):    
    # Initialize the jobscrparams dictionary with mandatory/useful keywords.
    jobscrparams={'memory': '4G',
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
    for job in jobsdir:
        jobs = job.rsplit('/',1)
        jobname,jobtrunc,jobfile=job_check(jobs[1],job)
        ld_extra,tc_exe,modname=tc_check(args.terachem)
        output=open(job+'/'+'jobscript','w')
        output.write('#$ -S /bin/bash\n')
        output.write('#$ -N %s\n' %(jobtrunc))
        output.write('#$ -l gpus=1\n')
        output.write('#$ -cwd\n')
        output.write('#$ -R y\n')
        output.write('#$ -l h_rt=%s\n' %(jobscrparams['time']))
        output.write('#$ -l h_rss=%s\n' %(jobscrparams['memory']))      
        output.write('#$ -q %s\n' %(jobscrparams['queue']))
        output.write('#$ -pe smp %s\n' %(jobscrparams['gpus']))
        output.write('# -fin terachem_input\n')
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
        output.write('terachem terachem_input > $SGE_O_WORKDIR/%s.out\n' %(jobname))
        output.close()
