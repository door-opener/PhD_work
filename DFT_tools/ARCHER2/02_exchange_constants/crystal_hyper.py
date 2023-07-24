import os
import sys
import numpy as np
from shutil import copy
import subprocess

def set_angle(cell_matrix, delta_angle, conv_matrix):
    # structure: primitive cell structure
    # delta_angle: delta on the conventional cell angle
    # conv_matrix: conversion matrix from the primitive to conventional cell

    #delta_angle = 1.
    #shift = 10.
    to_rad = np.pi/180
    to_deg = 180/np.pi

    conv_matrix_inv = np.linalg.inv(conv_matrix)
    #cubic_lattice = np.matmul(conv_matrix,structure.lattice.matrix)
    cubic_lattice = np.matmul( conv_matrix, cell_matrix )

    A = cubic_lattice[0,:]
    B = cubic_lattice[1,:]
    C = cubic_lattice[2,:]

    len_A = np.linalg.norm(cubic_lattice[0,:])
    len_B = np.linalg.norm(cubic_lattice[1,:])
    len_C = np.linalg.norm(cubic_lattice[2,:])

    A_norm = A/len_A
    B_norm = B/len_B
    C_norm = C/len_C

    alpha = np.arccos(np.dot(B_norm,C_norm))*(180/np.pi)
    beta = np.arccos(np.dot(A_norm,C_norm))*(180/np.pi)
    gamma = np.arccos(np.dot(A_norm,B_norm))*(180/np.pi)

    vector_sum = A_norm + B_norm + C_norm
    len_vector_sum = np.linalg.norm(vector_sum)

    vector_sum_norm = vector_sum/len_vector_sum

    D_a = np.cross(np.cross(A_norm,vector_sum_norm),A_norm)
    D_b = np.cross(np.cross(B_norm,vector_sum_norm),B_norm)
    D_c = np.cross(np.cross(C_norm,vector_sum_norm),C_norm)

    len_D_a = np.linalg.norm(D_a)
    len_D_b = np.linalg.norm(D_b)
    len_D_c = np.linalg.norm(D_c)

    D_a_norm = D_a/len_D_a
    D_b_norm = D_b/len_D_b
    D_c_norm = D_c/len_D_c

    alpha_t = np.arccos(np.dot(A_norm,vector_sum_norm))*(180/np.pi)
    beta_t = np.arccos(np.dot(B_norm,vector_sum_norm))*(180/np.pi)
    gamma_t = np.arccos(np.dot(C_norm,vector_sum_norm))*(180/np.pi)

    #len_vector_sum_p_tmp  = np.sqrt((3+6*(np.cos((alpha+delta_angle)*to_rad))))
    len_vector_sum_p = np.sqrt(3 + 2*(np.cos((alpha+delta_angle)*to_rad) +
                                        np.cos((beta+delta_angle)*to_rad) +
                                        np.cos((gamma+delta_angle)*to_rad)))


    alpha_p = (np.arccos((  1+ np.cos((beta+delta_angle)*to_rad) + np.cos((gamma+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg
    beta_p = (np.arccos((  1+ np.cos((alpha+delta_angle)*to_rad) + np.cos((gamma+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg
    gamma_p = (np.arccos((  1+ np.cos((alpha+delta_angle)*to_rad) + np.cos((beta+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg


    delta_alpha_t = (alpha_t-alpha_p)*to_rad
    delta_beta_t = (beta_t-beta_p)*to_rad
    delta_gamma_t = (gamma_t-gamma_p)*to_rad

    #print(alpha,alpha_t,alpha_p,delta_alpha_t)

    A1 = (np.cos(delta_alpha_t)*A_norm + np.sin(delta_alpha_t)*D_a_norm) #* len_A
    B1 = (np.cos(delta_beta_t)*B_norm + np.sin(delta_beta_t)*D_b_norm) #* len_B
    C1 = (np.cos(delta_gamma_t)*C_norm + np.sin(delta_gamma_t)*D_c_norm) #* len_C

    new_lattice_conv = np.array([A1,B1,C1])

    new_lattice_conv[0,:] = new_lattice_conv[0,:] *(len_A)
    new_lattice_conv[1,:] = new_lattice_conv[1,:] *(len_B)
    new_lattice_conv[2,:] = new_lattice_conv[2,:] *(len_C)
    return np.matmul(conv_matrix_inv,new_lattice_conv)

def get_angle(vect1, vect2):
    m1 = np.sqrt( np.power( vect1[0], 2.0) + np.power( vect1[1], 2.0) + np.power( vect1[2], 2.0) )
    m2 = np.sqrt( np.power( vect2[0], 2.0) + np.power( vect2[1], 2.0) + np.power( vect2[2], 2.0) )
    return np.arccos( ( ( np.dot(vect1, vect2) / (m1*m2) ) ) )

def getBasis(species, basis):
    basis_set = []
    repo = "/work/e05/e05/isa/05_CRYSTAL/00_UTILS/01_BASIS_SETS/"
    #repo = "/home/uccaitl/08_CRYSTAL/00_UTILS/02_GUESSDUAL_BASIS_SETS/"
    path = repo + [ x for x in  os.listdir( repo ) if species in x ][0]
    basis_file = path + "/" + [ x for x in os.listdir( path ) if basis in x ][0]
    f=open(basis_file, "r")
    [ next(f) for x in range(5) ]
    for line in f:
        split = line.split()
        basis_set.append( line )
    f.close()
    return basis_set

def writeInput(in_file, basis, pos, sym, spin, kpts, opt=False, CPHF=False, res=False):

    opt_args = [ "OPTGEOM", "FULLOPTG", "FRACTION", "PRINTFORCES", "PRINTOPT" ]
    dft_args = [ "DFT", "SPIN", "EXCHANGE", "PBESOL", "CORRELAT", "PBESOL", "HYBRID", "25", "TOLLDENS", "9" ]
    scf_btm = [ "FMIXING", "20", "SLOSHING", "SLOSHFAC", "1.5", "MAXCYCLE", "500", "TOLDEE", "10" ]

    if res:
       scf_top = [ "GUESSP", "TOLINTEG", "12 12 12 12 24", "SCFDIR"]
    else:
       scf_top = [ "TOLINTEG", "12 12 12 12 24", "SCFDIR"]

    CPHF_args = [ "CPHF", "FMIXING", "30" ]

    for i in range(len(pos)):
        in_file.write("%d  %12.14f  %12.14f  %12.14f\n" %(sym[i], pos[i][0], pos[i][1], pos[i][2]) )
    if opt:
       [ in_file.write("%s\n" %(x) ) for x in opt_args ]
       in_file.write("ENDOPT\n")
    if CPHF:
       [ in_file.write("%s\n" %(x) ) for x in opt_args ]
    in_file.write("ENDG\n")
    #
    for i in range(len(basis)):
        [ in_file.write(x) for x in basis[i] ]
    in_file.write("99 0\nENDBS\n")
    #
    [ in_file.write("%s\n"%(x)) for x in dft_args ]
    in_file.write("ENDDFT\n")
    #
    [ in_file.write("%s\n" %(x)) for x in scf_top ]
    [ in_file.write("%s\n" %(x)) for x in kpts ]
    [ in_file.write("%s\n" %(x)) for x in scf_btm ]
    #
    [ in_file.write("%s\n" %(x) ) for x in spin ]
    in_file.write("ENDSCF")

def SubmitJob(args, run=False, restart=False, short=False):
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
    if short:
       f.write("#SBATCH --time=00:20:00\n\n")
    else:
       f.write("#SBATCH --time=%d:00:00\n\n" %(int(args[1])) )
    f.write("#SBATCH --account=e05-gc-smw\n")
    f.write("#SBATCH --partition=standard\n")
    if short:
       f.write("#SBATCH --qos=short\n")
    else:
       f.write("#SBATCH --qos=standard\n")
    f.write("#SBATCH --export=none\n")
    f.write("\n")
    f.write("module load epcc-job-env\n")
    f.write("module load other-software\n")
    f.write("module load crystal\n\n")
    f.write("export FI_MR_CACHE_MAX_COUNT=0\n")
    f.write("\n")
    if restart:
       f.write("timeout %dm /work/e05/e05/isa/runCRYSTAL/Pcry_slurm %s %s &\nwait\n" %( (int(args[1])*60)-3 , args[3].replace(".d12",""), args[3].replace(".d12","") ))
    else:
       f.write("timeout %dm /work/e05/e05/isa/runCRYSTAL/Pcry_slurm %s &\nwait\n" %( (int(args[1])*60)-3 , args[3].replace(".d12","") ))
    f.write("/work/e05/e05/isa/runCRYSTAL/post_proc_slurm crys %s &\nwait\n" %(args[3].replace(".d12","") ))
    f.close()

    if run:
       err = subprocess.check_output("sbatch jobSubmit.slurm", stderr=subprocess.STDOUT, shell=True).decode('utf-8')
       print("%s"%(err) )

def getSymm(a_param, m_ordering ):
    """
    Pass in a0 get correct a.
    """
    a_cubic = np.float64(a_param*2.0)
    a_trigonal = np.around( np.float64( a_cubic*(np.sqrt(1.5)) ), decimals=8)
    a_prim = np.around( np.float64(a_cubic/(np.sqrt(2.0))), decimals=8)

    AFII_pos, AFII_sym = [ [0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75] ], [ 28, 8, 28, 8 ]
    FM_pos, FM_sym = [ [0.0, 0.0, 0.0], [0.5, 0.5, 0.5] ], [28, 8]
    AFI_pos, AFI_sym = [ [0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0], [0.0, 0.0, 0.5] ], [ 28, 28, 8, 8 ]
    AF3_pos, AF3_sym = [ [ 0.00000, 0.75000, 0.12500], [0.00000, 0.25000, 0.37500], [0.0000, 0.25000, 0.62500] ], [ 25, 25, 8 ]
    AF4_pos, AF4_sym = [ [ 0.25000, 0.25000, 0.12500], [0.25000, 0.25000, 0.62500], [0.2500, 0.75000, 0.37500], [0.25000, 0.75000, 0.87500] ], [25, 25, 8, 8]

    cell, setting, spg, spin, pos, sym, kpts = None, None, None, None, None, None, None

    if m_ordering == "AFII":
       setting = [0, 1, 0]
       spin = [ "ATOMSPIN", "2", "1 +1 3 -1", "SPINLOCK", "0 5" ]
       spg = 160
       pos, sym = AFII_pos, AFII_sym
       kpts = [ "SHRINK","17 34" ]
       cell = [ a_trigonal ]

    elif m_ordering == "AFI":
       setting = [0, 0, 0]
       spin = [ "ATOMSPIN", "2", "1 +1 2 -1", "SPINLOCK", "0 5" ]
       spg = 123
       cell = [ a_prim, a_cubic ]
       pos, sym = AFI_pos, AFI_sym
       kpts = [ "SHRINK", "0 38", "19 19 15" ]

    elif m_ordering == "FM":
       setting = [0, 1, 0]
       spin = [ "ATOMSPIN", "2", "1 +1 3 +1", "SPINLOCK", "4 5" ]
       spg = 160
       pos, sym = AFII_pos, AFII_sym
       kpts = [ "SHRINK","17 34" ]
       cell = [ a_trigonal ]

    elif m_ordering == "AF3":
       setting = [0, 0, 0]
       spin = [ "ATOMSPIN", "2", "1 +1", "1 -1", "SPINLOCK", "0 5" ]
       spg = 141
       cell = [ a_cubic, np.around(a_cubic*2.0, decimals=8) ]
       pos, sym = AF3_pos, AF3_sym
       kpts = [ "SHRINK", "0 26", "13 13 7" ]

    elif m_ordering == "AF4":
       setting = [0, 0, 0]
       spin = [ "ATOMSPIN", "2", "1 +1", "1 -1", "SPINLOCK", "0 5" ]
       spg = 59
       cell = [ a_cubic, a_prim, 2.0*a_prim ]
       pos, sym = AF4_pos, AF4_sym
       kpts = [ "SHRINK", "0 38", "13 19 9" ]

    return setting, spin, spg, cell, pos, sym, kpts

def main():

    Mn_tag, O_tag = "Ni_Kr.basis", "O_pob_TZVPP.basis"
    O_basis = getBasis("08_O", O_tag)
    Mn_basis = getBasis("28_Ni", Mn_tag)
    basis = [ Mn_basis, O_basis ]

    m_ordering = [ "FM", "AFI", "AFII" ]

    walltime, nprocs = 1, 16

    workdir = "01_PBESOL0/"

    os.mkdir(workdir)
    os.chdir(workdir)

    a_sta, a_end, numpts = 2.00, 2.15, 30
    da = (a_end - a_sta)/numpts

    for j in range(len(m_ordering)):

        outer_cur = ("%.2d_%s/"%(j+1, m_ordering[j]))
        os.mkdir(outer_cur)
        os.chdir(outer_cur)

        print(outer_cur)

        for i in range( numpts + 1 ):
            if i:
               a_cur = np.around(a_cur + da, decimals=3)
            else:
               a_cur = a_sta

            dir_cur = ("%.2d_%4.3fA/" %(i+1, a_cur) )
            job_cur = ("%s_%d"%((m_ordering[j], i+1)))
            print(" STEP %d " %(i+1) )
            job_args = [job_cur, walltime, nprocs, "NiO.d12" ]

            os.mkdir( dir_cur )
            os.chdir( dir_cur )

            if m_ordering[j] == "FM" or m_ordering[j] == "AFII":
               angle = 33.5573
            else:
               angle = None

            setting, spin, spg, cell, pos, sym, kpts = getSymm(a_cur, m_ordering[j])
            #setting, spin, spg, cell, pos, sym, angle = getSymm(2.10, "AFII", 90.0)

            #if m_ordering[j] == "AFII":
            #   restart = True
            #   short = True
            #else:

            short = False
            restart = False
            f=open("NiO.d12", "w")
            f.write("MnO\n")
            f.write("CRYSTAL\n")
            f.write("%d %d %d\n" %(setting[0], setting[1], setting[2] ))
            f.write("%d\n" %(spg) )
            if m_ordering[j] == "FM" or m_ordering[j] == "AFII": 
               f.write("%10.8f  %10.8f\n"%(cell[0], angle))
            elif m_ordering[j] == "AFI" or m_ordering[j] == "AF3":
               f.write("%10.8f  %10.8f\n"%(cell[0], cell[1]))
            else:
               f.write("%10.8f  %10.8f  %10.8f\n"%(cell[0], cell[1], cell[2]))
            f.write("%d\n" %(len(pos)))
            writeInput(f, basis, pos, sym, spin, kpts, opt=False, CPHF=False, res=False)
            #writeExt( cell, pos, sym )
            f.close()

            SubmitJob( job_args, run=False, restart=False, short=False)

            os.chdir("../")

        os.chdir("../")

if __name__=="__main__":
   main()
