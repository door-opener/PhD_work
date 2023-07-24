import os
import sys
import numpy as np
from shutil import copy
import subprocess

def SubmitJob(args, run=False):
    """
    args[0] = Job Name
    args[1] = walltime (hrs)
    args[2] = nprocs
    args[3] = input file name (CRYSTAL)
    """
    #Write jobscript.
    f=open("jobSubmit.slurm", "w")
    f.write("#!/bin/bash -l\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks-per-node=%d\n"%(args[2]))
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --job-name=%s\n" %(args[0]) )
    f.write("#SBATCH --time=%d:00:00\n\n" %(int(args[1])) )
    f.write("#SBATCH --account=e05-gc-smw\n")
    f.write("#SBATCH --partition=standard\n")
    f.write("#SBATCH --qos=standard\n")
    f.write("#SBATCH --export=none\n")
    f.write("\n")
    f.write("module load epcc-job-env\n")
    f.write("module load other-software\n")
    f.write("module load crystal\n\n")
    f.write("export FI_MR_CACHE_MAX_COUNT=0\n")
    f.write("\n")
    f.write("timeout %dm /work/e05/e05/isa/runCRYSTAL/Pcry_slurm %s &\nwait\n" %( (int(args[1])*60)-3 , args[3].replace(".d12","") ))
    f.write("/work/e05/e05/isa/runCRYSTAL/post_proc_slurm crys %s &\nwait\n" %(args[3].replace(".d12","") ))
    f.close()

    if run:
       err = subprocess.check_output("sbatch jobSubmit.slurm", stderr=subprocess.STDOUT, shell=True).decode('utf-8')
       print("%s"%(err) )

def main():

    try:
      in_dir = sys.argv[1]
    except:
      print("\nENTER DIR\n")
      sys.exit()

    dirs = sorted([ x for x in os.listdir( in_dir ) if "0" in x ] )

    for i in range(len( dirs )):

        os.chdir( dirs[i] )

        SubmitJob( [ ("DUMMY%d"%(i+1)), 1, 8, "NiO" ] )

        os.chdir( "../" )

if __name__ == "__main__":
   main()
