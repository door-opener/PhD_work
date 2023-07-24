import os, sys
import numpy as np
import subprocess as sp

def set_angle(cell_matrix, delta_angle, conv_matrix):
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

def getData( in_file ):
    fit_data = np.zeros( (47, 4), dtype=np.float64 )
    f=open(in_file, "r")
    for i,line in enumerate(f):
        split = line.split()
        fit_data[i][0] = split[0]
        fit_data[i][1] = split[1]
        fit_data[i][2] = split[2]
        fit_data[i][3] = split[3]
    f.close()
    return fit_data

def callGULP(infile):
    gulp_exe = "/home/uccaitl/06_Software/04_GULP_6/gulp-6.0_spingulp/Src/gulp"
    args = ("%s < %s > %s" %(gulp_exe, infile, infile.replace(".gin",".gout") ))
    err = sp.check_output(("%s" %(args)), stderr=sp.STDOUT, shell=True).decode('utf-8')
    return err

def writeGIN(keywords, option_words, pos, sym, cell, fname):
    potlib = "/home/uccaitl/07_GULP_SIMS/00_POTS/02_SET2/MnOc.lib"
    f=open(fname, "w")
    f.write("#\n#\n")
    [ f.write("%s " %(x) ) for x in keywords ]
    f.write("\n")
    #f.write("cell\n")
    #f.write("%6.5f  %6.5f  %6.5f  %6.4f  %6.4f  %6.4f\n" %(geom_in[0][0], geom_in[0][1], geom_in[0][2], geom_in[0][3], geom_in[0][4], geom_in[0][5]) )
    f.write("vectors\n")
    [ f.write("  %10.8f %10.8f %10.8f\n"%(x[0], x[1], x[2])) for x in cell ]
    f.write
    f.write("fractional\n")
    for i in range(len(pos)):
        f.write("%s %s %10.8f %10.8f %10.8f\n" %(sym[i], "core", pos[i][0], pos[i][1], pos[i][2]) )
    f.write("#\n")
    f.close()
    os.system("cat %s >> %s" %(potlib, fname) )
    f=open(fname, "a")
    f.write("#\n")
    f.write("maxcyc %d\n" %(1000) )
    f.write("xtol %2.1f\n" %(7.0) )
    f.write("gtol %2.1f\n" %(7.0) )
    f.write("ftol %2.1f\n" %(7.0) )
    [ f.write("%s\n" %(x)) for x in option_words ]
    f.close()

def writeGIN_fit( fit_data, pos, sym, fname ):

    potlib = "/home/uccaitl/07_GULP_SIMS/00_POTS/05_SET5/MnO_3.lib"
    keywords = [ "fit", "noflag", "comp" ]

    f=open(fname, "w")
    [ f.write("%s " %(x) ) for x in keywords ]
    f.write("\n")

    c_mat = np.asarray( [ [ -1.5, 0.5, 0.5 ], [ 0.5, -1.5, 0.5 ], [ 0.5, 0.5, -1.5 ] ], dtype=np.float64 )
    #s_a, s_angle = 1.0144, 0.9890
    s_a, s_angle = 1.0, 1.0

    for i in range(len( fit_data )):
        a_cur, alp_cur, E_cur, W_cur = fit_data[i][0], fit_data[i][1], fit_data[i][2], fit_data[i][3]
        a, alp = np.around( a_cur*s_a, decimals=8 ), np.around( alp_cur*s_angle, decimals=8 )
        cell_cur = getCell( 2.0*a )
        cell_cur = set_angle( cell_cur, (alp-90.0), c_mat )

        f.write("# Structure #%d | a = %8.6f A  alp = %8.6f deg.\n" %(i, a, alp) )
        f.write("vectors\n")
        [ f.write("  %10.8f %10.8f %10.8f\n"%(x[0], x[1], x[2])) for x in cell_cur ]
        f.write("fractional\n")
        for j in range(len(pos)):
            f.write("%s %s %10.8f %10.8f %10.8f\n" %(sym[j], "core", pos[j][0], pos[j][1], pos[j][2]) )
        f.write("#\n")
        f.write("observable\n")
        f.write("energy eV\n")
        f.write("%8.6f %f\n" %(E_cur, W_cur))
        f.write("end\n")
    f.close() 
    os.system("cat %s >> %s" %(potlib, fname) )

def getCell( a_param ):
    xform = np.asarray( [ [ 1.0, 0.5, 0.5 ], [ 0.5, 1.0, 0.5 ], [ 0.5, 0.5, 1.0 ] ], dtype=np.float64 )
    cubic = np.asarray( [ [ a_param, 0.0, 0.0, ], [ 0.0, a_param, 0.0 ], [ 0.0, 0.0, a_param ] ], dtype=np.float64 )
    return np.matmul( xform, cubic )

def main():

    try:
      in_file = sys.argv[1]
    except:
      print("\nENTER A FILE\n")
      sys.exit()

    fit_data = getData( in_file )
    pos = np.asarray( [ [ 0.0, 0.0, 0.0 ], [ 0.25, 0.25, 0.25 ], [ 0.5, 0.5, 0.5 ], [ 0.75, 0.75, 0.75 ] ], dtype=np.float64 )
    sym = [ "Mn1", "O", "Mn2", "O" ]

    writeGIN_fit( fit_data, pos, sym, "test_fit1_weighted.gin" )

    sys.exit()

    for i in range(numpts + 1):
        if i:
          cur = np.around( cur + da, decimals=3 )
        else:
          cur = st

        keywords = [ "fit", "noflag" ]
        option_words = [ ]
        cell_cur = getCell( 2.0*cur )
        pos = np.asarray( [ [ 0.0, 0.0, 0.0 ], [ 0.25, 0.25, 0.25 ], [ 0.5, 0.5, 0.5 ], [ 0.75, 0.75, 0.75 ] ], dtype=np.float64 )
        sym = [ "Mn", "O", "Mn", "O" ]
        writeGIN( keywords, option_words, pos, sym, cell_cur, ("step%.2d.gin"%(i+1)) )

    sys.exit()

    callGULP( "dummy.gin" )

if __name__ == "__main__":
   main()
