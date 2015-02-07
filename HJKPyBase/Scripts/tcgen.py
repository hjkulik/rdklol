''' 
Created on 12/27/14

@author: hjk
'''
import argparse
import glob
import sys
import os
import subprocess
execfile('/usr/share/Modules/init/python.py')

def mybash(cmd):
    """ Tim's approach to getting the stdout from a bash command.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        stdout.append(line)        
        if line == '' and p.poll() != None:
            break
    return ''.join(stdout)

def tc_check(userinput):
    """ This checks for whether the user requested a good terachem install. If not it quits.
        It also returns the LD_LIBRARY_PATH spot if needed (for post-intbox TC) and the executable.
    """
    modname=''
    if not userinput:
      module('load','terachem/alpha')
      modname='terachem/alpha'
    else:
      if userinput[0] is not '/':
         if mybash('/usr/bin/modulecmd python avail %s' %(userinput)):
            module('unload','terachem')
            module('load',userinput)
            modname=userinput
         elif mybash('/usr/bin/modulecmd python avail terachem/%s' %(userinput)):
            module('unload', 'terachem')
            modname='terachem/'+userinput
            module('load','terachem/'+userinput)
         else:
            print("You're pointing to a modulefile that doesn't exist or not specifying a full path: %s. Try again!" %(userinput))
            exit()
      else:
         os.environ['TeraChem'] = userinput
      if not glob.glob(os.path.expandvars('$TeraChem/basis/')):
         print("Custom TeraChem directory %s does not have basis files. Try again!" %(userinput))
         exit()
    ld_extra=''
    if glob.glob(os.path.expandvars('$TeraChem/lib/')):
       ld_extra=glob.glob(os.path.expandvars('$TeraChem/lib/'))[0]
    if glob.glob(os.path.expandvars('$TeraChem/bin/terachem')):
       tc_exe=glob.glob(os.path.expandvars('$TeraChem/bin/terachem'))[0]
    elif glob.glob(os.path.expandvars('$TeraChem/bin/int')):
       tc_exe=glob.glob(os.path.expandvars('$TeraChem/bin/int'))[0]
    elif glob.glob(os.path.expandvars('$TeraChem/int')):
       tc_exe=glob.glob(os.path.expandvars('$TeraChem/int'))[0]
    elif glob.glob(os.path.expandvars('$TeraChem/terachem')):
       tc_exe=glob.glob(os.path.expandvars('$TeraChem/terachem'))[0]
    else:
       print("TeraChem directory %s does not have a correctly named terachem executable. Try again!" %(glob.glob(os.path.expandvars('$TeraChem'))))
       exit()
    return ld_extra,tc_exe,modname

def job_check(userinput):
    """ Checks for existence of a pdb/xyz file based on job file existing.
        Returns truncated jobname (e.g. for logging/queue script) and jobfile name.
    """
    if userinput:
       jobname=userinput
       if glob.glob(jobname+'.xyz'):
          jobfile=jobname+'.xyz'
       elif glob.glob(jobname+'.pdb'):
          jobfile=jobname+'.pdb'
       else:
          print("No pdb/xyz file with prefix %s found! Try again." %(userinput))
          exit()
       if len(jobname) > 8:
          if '_' in jobname:
             jobtrunc=userinput.split('_')[0][0:4]+userinput.split('_')[1][0:4]
          else:
             jobtrunc=userinput[0:8]
       else:
          jobtrunc=jobname
    else:
       print("You didn't specify a jobname. This is needed so we know which xyz/pdb to use. Try again.")
       exit()
    return jobname,jobtrunc,jobfile

if __name__ == "__main__":
    # Parse commandline arguments
    parser= argparse.ArgumentParser()
    parser.add_argument("-j","--job", help="prefix for pdb/xyz structure file and jobname.")
    parser.add_argument("-m","--method", help="electronic structure approach for terachem job. Specify UHF/UDFT to enforce levelshifting and unrestricted calculation for OS singlets (default: b3lyp).")
    parser.add_argument("-b","--basis", help="basis for terachem job (default: LACVP*).") 
    parser.add_argument("-s","--spin", help="spin for system (default: singlet).")
    parser.add_argument("-ch","--charge", help="charge for system (default: neutral).") 
    parser.add_argument("-x","--extra", nargs="+", help="extra arguments in syntax tckeyword=value. suggestions: nstep (opt), maxit (scf), scf=diis+a, min_tolerance>4.5e-4/min_tolerance_e>1.0e-6, min_coordinates=cartesian, orbitalswrtfrq=<num>,multibasis=<basis file>.")
    parser.add_argument("-tc","--terachem", help="full path to custom terachem installation.")
    # TBD, could do groupings of smart keywords at some point - eg 'floppy' or 'mustconverge' but not for now.
    args=parser.parse_args()
    parser.parse_args()
    # Initialize the jobparams dictionary with mandatory/useful keywords.
    jobparams={'run': 'minimize',
           'timings': 'yes',
           'nstep': '1001', 
           'min_method': 'bfgs', 
           'min_init_hess': 'fischer-almlof',
           'min_tolerance': '4.5e-4',
           'min_tolerance_e': '1.0e-6',
           'maxit': '150',
           'scrdir': './scr',
           'method': 'b3lyp',
           'basis': 'lacvps_ecp',
           'spinmult': '1',
           'charge': '0'}
    # Overwrite plus add any new dictionary keys from commandline input.       
    # Check to make sure that $TeraChem is loaded/initialized:
    ld_extra,tc_exe,modname=tc_check(args.terachem)
    # Setting jobname for files + truncated name for queue.
    jobname,jobtrunc,jobparams['coordinates']=job_check(args.job)
    # Method parsing, does not check if a garbage method is used here:
    unrestricted=False
    if args.method:
       if ('u' or 'U') in args.method[0]:
           #Unrestricted calculation
           jobparams['method']=args.method 
           unrestricted=True
       else:
           #Restricted calculation
           unrestricted=False
           jobparams['method']=args.method
           if args.spin and int(args.spin) > 1:
              jobparams['method']='u'+args.method
              unrestricted=True
    # Just carry over spin and charge keywords if they're set. Could do checks, none for now.
    if args.spin:
       jobparams['spinmult']=args.spin
    if args.charge:
       jobparams['charge']=args.charge
    # Check for existence of basis and sanitize name
    if args.basis:
       ecp=False # Flag not currently used, for deciding gpus_ecp code or not later. Can always specify with 'extra' command
       if '*' in args.basis:
          jobparams['basis']=args.basis.replace('*','s')
       else:
          jobparams['basis']=args.basis
       if 'ecp' in args.basis:
          ecp=True
       if not glob.glob('%s/basis/%s' %(os.path.expandvars('$TeraChem'),basis)): 
          print("Basis %s does not exist. Try again!" %(args.basis))
          exit()
    # Overwrite plus add any new dictionary keys from commandline input.       
    if args.extra:
       for elem in range(0,len(args.extra)):
          key,val=args.extra[elem].split('=')
          jobparams[key]=val
    # Extra keywords for unrestricted. 
    if unrestricted:
       # If running unrestricted, assume convergence will be more difficult for now.
       jobparams['scf']='diis+a' 
       if not jobparams.has_key('levelshift'):
          jobparams['levelshift']='yes'
       elif jobparams['levelshift'] != 'yes':
          print("Warning! You're doing an unrestricted calculation but have set levelshift = %s" %(jobparams['levelshift']))
       if not jobparams.has_key('levelshiftvala'):
          jobparams['levelshiftvala']='1.6'
       if not jobparams.has_key('levelshiftvalb'):
          jobparams['levelshiftvalb']='0.1'
    # Now we're ready to start building the input file and the job script
    output=open(jobname+'.in','w')
    for keys in jobparams.keys():
        output.write('%s %s\n' %(keys,jobparams[keys]))
    output.write('end\n')
    output.close()
