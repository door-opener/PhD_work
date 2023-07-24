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

def writeInput(in_file, basis, pos, sym, spin, kpts, opt=False, CPHF=False):

    opt_args = [ "OPTGEOM", "FULLOPTG", "FRACTION", "PRINTFORCES", "PRINTOPT" ]
    dft_args = [ "DFT", "SPIN", "EXCHANGE", "PBESOL", "CORRELAT", "PBESOL", "HYBRID", "25", "TOLLDENS", "9" ]
    scf_top = [ "TOLINTEG", "12 12 12 12 24", "SCFDIR"]
    scf_btm = [ "FMIXING", "20", "SLOSHING", "SLOSHFAC", "1.5",  "MAXCYCLE", "500", "TOLDEE", "10" ]
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

def writeExt( cell, pos, sym ):
    f=open("MnO.gui", "w")
    f.write("3 1 1\n")
    [ f.write("  %10.8f %10.8f %10.8f\n" %(x[0], x[1], x[2]) ) for x in cell ]
    f.write("1\n")
    [ f.write("  %3.2f %3.2f %3.2f\n" %(x[0], x[1], x[2]) ) for x in np.identity(3) ]
    f.write("  0.00 0.00 0.00\n")
    f.write("%d\n" %(len(pos)) )
    for i in range(len(pos)):
        f.write("  %d  %14.12f  %14.12f  %14.12f\n" %(sym[i], pos[i][0], pos[i][1], pos[i][2]))
    f.close()

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

def getSymm(a_param, m_ordering, dalp):
    """
    Pass in a0 get correct a.
    """
    a_cubic = np.float64(a_param*2.0)
    a_trigonal = np.around( np.float64( a_cubic*(np.sqrt(1.5)) ), decimals=8)
    a_prim = np.around( np.float64(a_cubic/(np.sqrt(2.0))), decimals=8)
    conv_mat = np.asarray( [ [ -1.5, 0.5, 0.5 ], [ 0.5, -1.5, 0.5 ], [ 0.5, 0.5, -1.5 ] ], dtype=np.float64 )

    AFII_pos, AFII_sym = [ [0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75] ], [ 25, 8, 25, 8 ]
    cell, setting, spg, spin, pos, sym = None, None, None, None, None, None

    if m_ordering == "AFII":
       setting = [0, 1, 0]
       spin = [ "ATOMSPIN", "2", "1 +1 3 +1", "SPINLOCK", "10 8" ]
       spg = 160
       pos, sym = AFII_pos, AFII_sym

       # Generate standard primitive #
       dummy = np.asarray( [ [ a_cubic, a_param, a_param ], [ a_param, a_cubic, a_param ], [ a_param, a_param, a_cubic ] ], dtype=np.float64 )

       # Adjust angle #
       cell = set_angle( dummy, dalp, conv_mat )
       angle = get_angle(cell[0], cell[1])*180.0/np.pi

    return setting, spin, spg, cell, pos, sym, angle

def main():

    Mn_tag, O_tag = "Mn_Karlsruhe_def_TZVP.basis", "O_pob_TZVPP.basis"
    O_basis = getBasis("08_O", O_tag)
    Mn_basis = getBasis("25_Mn", Mn_tag)
    basis = [ Mn_basis, O_basis ]

    walltime, nprocs = 2, 16

    workdir = "PBESOL0_GRID_FM/"

    os.mkdir(workdir)
    os.chdir(workdir)

    a_sta, a_end, numpts = 2.10, 2.40, 10
    alp_st, alp_end = 90.0, 91.5

    da = (a_end - a_sta)/numpts
    dalp = (alp_end - alp_st)/numpts

    for i in range( numpts + 1 ):
        if i:
           a_cur = np.around(a_cur + da, decimals=3)
        else:
           a_cur = a_sta

        for j in range( numpts + 1 ):
            if j:
               alp_cur = np.around(alp_cur + dalp, decimals=3 )
            else:
               alp_cur = alp_st

            dir_cur = ("%d_%3.2fA_%3.2fALP/" %((11*i+j)+1, a_cur, alp_cur) )
            job_cur = ("AFII_%d"%((11*i+j)+1))
            print(" STEP %d " %(11*i+j) )
            job_args = [job_cur, walltime, nprocs, "MnO.d12" ]

            os.mkdir( dir_cur )
            os.chdir( dir_cur )

            setting, spin, spg, cell, pos, sym, angle = getSymm(a_cur, "AFII", np.around(alp_cur - 90.0, decimals=4) )
            #setting, spin, spg, cell, pos, sym, angle = getSymm(2.10, "AFII", 90.0)
            kpts = [ "SHRINK","17 34" ]
            a = np.around( np.linalg.norm( cell[0] ), decimals=12 )
            angle = np.around( angle, decimals=12 )

            f=open("MnO.d12", "w")
            f.write("MnO\n")
            f.write("CRYSTAL\n")
            f.write("%d %d %d\n" %(setting[0], setting[1], setting[2] ))
            f.write("%d\n" %(spg) )
            f.write("%14.12f  %14.12f\n"%(a, angle))
            f.write("%d\n" %(len(pos)))
            writeInput(f, basis, pos, sym, spin, kpts, opt=False)
            #writeExt( cell, pos, sym )
            f.close()
            SubmitJob( job_args, run=False )

            os.chdir("../")

if __name__=="__main__":
   main()
