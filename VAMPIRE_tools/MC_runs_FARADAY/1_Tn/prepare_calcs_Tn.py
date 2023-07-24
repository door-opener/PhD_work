import os, sys
import subprocess as sp
import numpy as np
from shutil import copy

def writeJob(name, run=False):
    f=open("jobSubmit.sh", "w")
    f.write("#!/bin/bash -l\n")
    f.write("#$ -l h_rt=48:00:00\n")
    f.write("#$ -l mem=5G\n")
    f.write("#$ -N %s\n"%(name))
    f.write("#$ -wd %s\n\n" %(os.getcwd()))
    f.write("/home/uccaitl/06_Software/vampire_build_intel2018/vampire-serial")
    f.close()
    if run:
       err = sp.check_output( "qsub jobSubmit.sh", stderr=sp.STDOUT, shell=True).decode('utf-8')
       return err

def getMat( xc , corr=False ):
    repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/01_MnO_AND_NiO/00_RUN_FILES/"

    if corr:
       repo = ("%s04_MAT_CORR/"%(repo))
    else:
       repo = ("%s05_MAT_UNCORR/"%(repo))

    if xc == "PBE":
       dummy = [ x for x in os.listdir( repo ) if xc in x and "PBE0" not in x and "PBESol" not in x ]
    elif xc == "PBESol0":
       dummy = [ x for x in os.listdir( repo ) if xc in x and "NiO" not in x and "CRYSTAL17" not in x ]
    else:
       dummy = [ x for x in os.listdir( repo ) if xc in x ]
    return ("%s%s"%(repo,dummy[0]))

def getUCF( J_type, xc, hyper=False, corr=False, J1_split=False ):
    repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/01_MnO_AND_NiO/00_RUN_FILES/"
    if corr:
       repo = ("%s01_UCF_J2_CORR/"%(repo))
    elif corr == False and hyper == True and J1_split == False:
       repo = ("%s03_UCF_J3_J4/"%(repo))
    elif corr == False and hyper == False and J1_split == False:
       repo = ("%s02_UCF_J2_UNCORR/"%(repo))
    elif corr == False and hyper == False and J1_split == True:
       repo = ("%s08_UCF_SPLIT_J1/"%(repo))

    if xc == "PBE":
       dummy = [ x for x in os.listdir( repo ) if J_type in x and xc in x and "PBESol" not in x and "PBE0" not in x ]
    elif xc == "PBESol":
       dummy = [ x for x in os.listdir( repo ) if J_type in x and xc in x and "PBESol0" not in x and "CRYSTAL17" not in x ]
    elif xc == "PBESol0":
       dummy = [ x for x in os.listdir( repo ) if J_type in x and xc in x and "CRYSTAL17" not in x ]
    else:
       dummy = [ x for x in os.listdir( repo ) if J_type in x and xc in x ]
    return ("%s%s"%(repo, dummy[0]) )

def writeMCInputSTATIC(matfile, ucf_file, size, Tinc, Tmax):
    Tmin = 0.0
    numpts = Tmax/Tinc
    
    eqm_steps, loop_steps = 5000, 20000
    #print("\nTmax/min = %d/%dK | EQM steps = %d | LOOP steps = %d | #pts = %d\n" %(Tmax, Tmin, eqm_steps, loop_steps, numpts) )
    output_args = ["output:temperature", "output:material-magnetisation", "output:mean-susceptibility", "output:mean-total-energy" ]
    config_args = [ "config:atoms", ("config:atoms-output-rate=%d"%((Tmax/Tinc)-1)) ]
    pbc = [ "create:periodic-boundaries-x", "create:periodic-boundaries-y", "create:periodic-boundaries-z" ]

    cool_fn = "linear"
    precond_steps = 1000
    increment = 5

    f=open("input", "w")
    f.write("dimensions:system-size-x=%d !nm\ndimensions:system-size-y=%d !nm\ndimensions:system-size-z=%d !nm\n\n" %(size[0], size[1], size[2]))
    [ f.write("%s\n"%(x)) for x in pbc ]
    f.write("\nmaterial:file=%s\n" %(matfile) )
    f.write('material:unit-cell-file="%s"\n'%(ucf_file))
    f.write("\nsim:temperature-increment=%6.5f\n" %(Tinc) )
    f.write("sim:time-steps-increment=%d\n" %(increment) )
    f.write("sim:cooling-function=%s\n" %(cool_fn) )
    f.write("sim:maximum-temperature=%8.6f\n" %(Tmax) )
    f.write("sim:minimum-temperature=%8.6f\n" %(Tmin) )
    f.write("sim:preconditioning-steps=%d\n" %(precond_steps) )
    f.write("sim:equilibration-temperature=%3.2f\n" %(Tmax) )
    f.write("sim:equilibration-time-steps=%d\n" %(eqm_steps) )
    f.write("sim:loop-time-steps=%d\n" %(loop_steps) )
    f.write("\nsim:program=%s\n" %("curie-temperature"))
    f.write("sim:integrator=%s\n\n" %("monte-carlo") )
    [ f.write("%s\n" %(x) ) for x in output_args ]
    f.write("\n")
    [ f.write("%s\n" %(x) ) for x in config_args ]
    f.close()

def runJobs(run_typ, size, functionals, corr, J1_split):
    hyper,cry17=None, None
    Tinc, Tmax = 0.0, 0.0
    for i in range(len( functionals )):
        dir_cur = ("%.2d_%s/" %(i+1, functionals[i]))
        Jname = ("%s_%s" %(run_typ, functionals[i]))

        print(Jname)

        if run_typ != "J1" and "CRYSTAL17" in functionals[i]:
           print("\nSkipping...")
           continue

        os.mkdir( dir_cur )
        os.chdir( dir_cur )

        if "CRYSTAL17" not in functionals[i] and run_typ != "J1":
            hyper=True
        else:
            hyper=False

        if "NiO" in functionals[i]:
           Tinc, Tmax = 20.0, 1000.0
        else:
           Tinc, Tmax = 5.0, 200.0

        UCF_cur = getUCF( run_typ, functionals[i], hyper=hyper, corr=corr, J1_split=J1_split ) 
        mat_cur = getMat( functionals[i], corr)

        dest1 = ("%s/%s" %(os.getcwd(), UCF_cur.split("/")[-1]))
        dest2 = ("%s/%s" %(os.getcwd(), mat_cur.split("/")[-1]))

        copy(UCF_cur, dest1)
        copy(mat_cur, dest2)

        writeMCInputSTATIC( mat_cur.split("/")[-1], UCF_cur.split("/")[-1], size, Tinc, Tmax )
        writeJob(Jname, run=False)

        os.chdir( "../" )

def main():

    try:
      in_dir = ("%s/%s"%(os.getcwd(), sys.argv[1]))
    except:
      print("\nENTER A DIRECTORY\n")
      sys.exit()

    # Uncorrected J2 #
    #run_typ, size = [ "J1" ], [15, 15, 15]
    #functionals = [ "PBE", "PBESol", "SCAN", "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17", "PBE0_NiO_CRYSTAL17", "PBESol0_NiO_CRYSTAL17" ]
    #functionals = [ "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17", "PBE0_NiO_CRYSTAL17", "PBESol0_NiO_CRYSTAL17" ]
    #corr = False
    #J1_split = False

    # Corrected J2 #
    run_typ, size = [ "J1" ], [15, 15, 15]
    functionals = [ "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17" ]
    corr = True
    J1_split = False

    # Uncorrected J2 - J1 split #
    #run_typ, size = [ "J1" ], [15, 15, 15]
    #functionals = [ "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17", "PBE0_NiO_CRYSTAL17", "PBESol0_NiO_CRYSTAL17" ]
    #corr = False
    #J1_split = False

    # J1, J2, J3, J4 #
    #run_typ, size = [ "J3", "J4" ], [15, 15, 15]
    #functionals = [ "PBE", "PBESol", "SCAN", "PBE0", "PBESol0" ]
    #corr = False

    print("\nPreparing Calcs in: [%s]\n"%(in_dir))
    print("J2 Correction = %s \n" %(corr))
    print("J1 Split = %s \n" %(J1_split))
    os.chdir( in_dir )

    for i in range(len( run_typ )):

        os.mkdir("%.2d_%s/" %(i+1, run_typ[i]) )
        os.chdir("%.2d_%s/" %(i+1, run_typ[i]) )
        print("RUN TYPE = %s " %(run_typ[i]) )
        runJobs( run_typ[i], size, functionals, corr, J1_split )

        os.chdir("../")

if __name__ == "__main__":
    main()
