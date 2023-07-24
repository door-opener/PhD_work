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

def getMat( xc , J_type, corr=False, J1234=False ):
    repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/"

    if corr:
       repo = ("%s04_MAT_J2_CORR/"%(repo))
    elif corr==False and "EXPT" in xc or corr == False and "PEPY" in xc:
       repo = ("%s06_MAT_EXPT/"%(repo))
    elif J1234== True and corr == False:
       repo = ("%s09_J1234_mat/"%(repo))

    if xc == "PEPY_1" or xc == "PEPY_2":
       dummy = [ x for x in os.listdir( repo ) if "PEPY_1974_1" in x ]
    elif xc == "PEPY_3" or xc == "PEPY_4":
       dummy = [ x for x in os.listdir( repo ) if "PEPY_1974_2" in x ]
    elif "LINES" in xc:
       dummy = [ x for x in os.listdir( repo ) if "LINES" in x ]
    elif J1234==True and corr==False:
       if xc == "PBE":
          dummy = [ x for x in os.listdir( repo ) if xc in x and "PBE0" not in x and "PBESol" not in x ]
       elif xc == "PBESol0":
          dummy = [ x for x in os.listdir( repo ) if xc in x and "NiO" not in x and "CRYSTAL17" not in x ]
       elif xc == "PBESol":
          dummy = [ x for x in os.listdir( repo ) if xc in x and "PBESol0" not in x ]
       else:
          dummy = [ x for x in os.listdir( repo ) if xc in x ]
    else:
       dummy = [ x for x in os.listdir( repo ) if xc.split("EXPT_")[-1] in x ]

    return ("%s%s"%(repo,dummy[0]))

def getUCF( J_type, xc, hyper=False, corr=False, J1_split=False, J1234=False ):

    if "NiO" in xc:
       repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/02_UCF_NiO/"
       dummy = [ x for x in os.listdir( repo ) if xc in x ]
    elif "PEPY" in xc:
       repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/07_UCFS_PEPY/"
       dummy = [ x for x in os.listdir( repo ) if xc in x ]
    elif "EXPT" in xc:
       repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/05_UCF_EXPT/"
       dummy = [ x for x in os.listdir( repo ) if xc.split("EXPT_")[-1] in x ]
    elif J1234:
       repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/08_J1234_UCF/"
       dummy = [ x for x in os.listdir( repo ) if xc in x and J_type in x ]
    else:
       repo = "/home/uccaitl/05_VAMPIRE_SIMS/10_FINAL_MC_SIMS/02_MnO_NiO_CW/00_RUN_FILES/01_UCF_J1/"
       dummy = [ x for x in os.listdir( repo ) if ".ucf" in x ]
    return ("%s%s"%(repo, dummy[0]) )

def writeMCInputSTATIC(Tcur, matfile, ucf_file, size):

    dT = 0.0001
    Tmax = Tcur + dT
    numpts = 200
    eqm_steps, loop_steps = 5000, 20000
    Tinc = dT/numpts

    output_args = ["output:temperature", "output:material-magnetisation", "output:mean-susceptibility", "output:mean-total-energy" ]
    config_args = [ "config:atoms", ("config:atoms-output-rate=%d"%((numpts)-1)) ]
    pbc = [ "create:periodic-boundaries-x", "create:periodic-boundaries-y", "create:periodic-boundaries-z" ]

    cool_fn = "linear"
    precond_steps = 1000
    increment = 5

    f=open("input", "w")
    f.write("dimensions:system-size-x=%d !nm\ndimensions:system-size-y=%d !nm\ndimensions:system-size-z=%d !nm\n\n" %(size[0], size[1], size[2]))
    [ f.write("%s\n"%(x)) for x in pbc ]
    f.write("\nmaterial:file=%s\n" %(matfile) )
    f.write('material:unit-cell-file="%s"\n'%(ucf_file))
    f.write("\nsim:temperature-increment=%10.8f\n" %(Tinc) )
    f.write("sim:time-steps-increment=%d\n" %(increment) )
    f.write("sim:cooling-function=%s\n" %(cool_fn) )
    f.write("sim:maximum-temperature=%8.6f\n" %(Tmax) )
    f.write("sim:minimum-temperature=%8.6f\n" %(Tcur) )
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

def runJobsCW(run_typ, size, functionals, corr, J1_split, J1234):
    hyper,cry17,Tscale=None, None, None
    Tinc, Tmax = 0.0, 0.0
    if run_typ == "J3":
       #Tn_dict = {"PBE":121.5,"PBESol":108.67,"SCAN":89.5837,"PBE0":48.7,"PBESol0":52.409}
       Tn_dict = {"PBE":126.09,"PBESol":114.95,"SCAN":92.04,"PBE0":49.8388,"PBESol0":53.491}
    elif run_typ == "J4":
       Tn_dict = {"PBE":104.329,"PBESol":100.681,"SCAN":75.8225,"PBE0":49.8388,"PBESol0":44.098}
    else:
       # J2 UNCORRECTED #
       Tn_dict = {"PBE":121.5,"PBESol":108.67,"SCAN":89.5837,"PBE0":48.7,"PBESol0":52.409,"PBE0_CRYSTAL17":58.81,"PBESol0_CRYSTAL17":67.96,"PBE0_NiO_CRYSTAL17":320,"PBESol0_NiO_CRYSTAL17":395,"EXPT_LINES_1":37.5, "EXPT_LINES_2":116, "EXPT_HUTCHINGS":408, "EXPT_SHANKER":375.5, "PEPY_1":36.2, "PEPY_2":131.5, "PEPY_3":45, "PEPY_4":115}

    Tn_scale = [1.5, 2, 2.5, 3, 3.5 ]

    for i in range(len( functionals )):

        dir_cur = ("%.2d_%s/" %(i+1, functionals[i]))

        Jname = ("%s_%s" %(run_typ, functionals[i]))

        print(Jname)

        if run_typ != "J1" and "CRYSTAL17" in functionals[i]:
           print("\nSkipping...")
           continue

        if "CRYSTAL17" not in functionals[i] and run_typ != "J1":
            hyper=True
        else:
            hyper=False

        UCF_cur = getUCF( run_typ, functionals[i], hyper=hyper, corr=corr, J1_split=J1_split, J1234=J1234 ) 
        mat_cur = getMat( functionals[i], run_typ, corr=corr, J1234=J1234)

        Tn_cur = [ v for k,v in Tn_dict.items() if k == functionals[i] ][0]

        #if "NiO" in functionals[i]:
        #    T_scale = 1.0
        #elif "HUTCHINGS" or "SHANKER" in functionals[i]:
        #    T_scale = 1.0
        #elif "PEPY" or "LINES" in functionals[i]:
        #    T_scale = 2.5
        #else:
        #    T_scale = 2.5

        #if "HUTCHINGS" in functionals[i]:
        #   T_scale = 1.0
        #elif "SHANKER" in functionals[i]:
        #   T_scale = 1.0
        #else:
        #   T_scale = 2.5

        T_scale = 2.5

        os.mkdir( dir_cur )
        os.chdir( dir_cur )

        for j in range(len(Tn_scale)):

            T_cur = (Tn_cur*Tn_scale[j])/T_scale
            inner_cur = ("%.2d_%.1f/" %(j+1, Tn_cur*Tn_scale[j]))

            os.mkdir( inner_cur )
            os.chdir( inner_cur )

            dest1 = ("%s/%s" %(os.getcwd(), UCF_cur.split("/")[-1]))
            dest2 = ("%s/%s" %(os.getcwd(), mat_cur.split("/")[-1]))

            copy(UCF_cur, dest1)
            copy(mat_cur, dest2)

            writeMCInputSTATIC(T_cur, mat_cur.split("/")[-1], UCF_cur.split("/")[-1], size)
            writeJob(Jname, run=True)

            os.chdir( "../" )

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
    ##functionals = [ "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17", "PBE0_NiO_CRYSTAL17", "PBESol0_NiO_CRYSTAL17" ]
    #corr = False
    #J1_split = False

    # EXPT #
    #run_typ, size = [ "J1" ], [15, 15, 15]
    #functionals = [ "EXPT_LINES_1", "EXPT_LINES_2", "EXPT_HUTCHINGS", "EXPT_SHANKER", "PEPY_1", "PEPY_2", "PEPY_3", "PEPY_4" ]
    #corr = False
    #J1_split = False

    # J1234 #
    run_typ, size = [ "J3", "J4" ], [15, 15, 15]
    functionals = [ "PBE", "PBESol", "SCAN", "PBE0", "PBESol0" ]
    #functionals = [ "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17", "PBE0_NiO_CRYSTAL17", "PBESol0_NiO_CRYSTAL17" ]
    corr = False
    J1_split = False
    J1234 = True

    # Corrected J2 #
    #run_typ, size = [ "J1" ], [15, 15, 15]
    #functionals = [ "PBE0", "PBESol0", "PBE0_CRYSTAL17", "PBESol0_CRYSTAL17" ]
    #corr = True
    #J1_split = False

    #TEST#
    #run_typ, size = [ "J1" ], [15, 15, 15]
    #functionals = [ "PBE0" ]
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
        runJobsCW( run_typ[i], size, functionals, corr, J1_split, J1234 )

        os.chdir("../")

if __name__ == "__main__":
    main()
