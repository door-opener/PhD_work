from CRYSTALpytools.crystal_io import *
from CRYSTALpytools.convert import *
import numpy as np
import sys, os

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import *

from os.path import join
import shutil as sh

from ase.visualize import view
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.periodic_table import *

from mendeleev import element

def get_delta_angle(structure, range_angle, conv_matrix=[]):
    # Return a list of delta angles from a range of %

    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    lattice = np.matmul(conv_matrix,structure.lattice.matrix)

    angle = Structure(lattice,[],[]).lattice.alpha

    delta_angle = []
    for i in np.arange(range_angle[0],range_angle[1]+range_angle[2]/2,range_angle[2]):
        delta_angle.append(angle*i)

    return np.round(delta_angle,8)

def get_delta_lattice(structure, range_lattice, conv_matrix=[]):
    # Return a list of delta lattice from a range of %

    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    lattice = np.matmul(conv_matrix,structure.lattice.matrix)

    lattice = Structure(lattice,[],[]).lattice.a

    delta_lattice = []
    for i in np.arange(range_lattice[0],range_lattice[1]+range_lattice[2]/2,range_lattice[2]):
        delta_lattice.append(lattice*i)

    return np.round(delta_lattice,8)

def set_angle(structure, delta_angle, conv_matrix=[]):
    # structure: primitive cell structure
    # delta_angle: delta on the conventional cell angle
    # conv_matrix: conversion matrix from the primitive to conventional cell
    
    import numpy as np
    from pymatgen.core.structure import Structure

    to_rad = np.pi/180
    to_deg = 180/np.pi


    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    conv_matrix_inv = np.linalg.inv(conv_matrix)
    lattice = np.matmul(conv_matrix,structure.lattice.matrix)


    A = lattice[0,:]
    B = lattice[1,:]
    C = lattice[2,:]

    len_A = np.linalg.norm(lattice[0,:])
    len_B = np.linalg.norm(lattice[1,:])
    len_C = np.linalg.norm(lattice[2,:])

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

    new_lattice_prim = np.matmul(conv_matrix_inv,new_lattice_conv)


    return Structure(new_lattice_prim,structure.atomic_numbers,structure.frac_coords)  

def set_cell(structure, delta_cell, conv_matrix=[]):
    # structure: primitive cell structure
    # delta_cell: delta on the conventional cell vectors
    # conv_matrix: conversion matrix from the primitive to conventional cell
    
    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    conv_matrix_inv = np.linalg.inv(conv_matrix)
    lattice = np.matmul(conv_matrix,structure.lattice.matrix)


    A = lattice[0,:]
    B = lattice[1,:]
    C = lattice[2,:]

    len_A = np.linalg.norm(lattice[0,:])
    len_B = np.linalg.norm(lattice[1,:])
    len_C = np.linalg.norm(lattice[2,:])


    ###
    lattice_new = np.zeros((3,3))
    lattice_new[0,:] = (lattice[0,:]/len_A) * (len_A+delta_cell)
    lattice_new[1,:] = (lattice[1,:]/len_B) * (len_B+delta_cell)
    lattice_new[2,:] = (lattice[2,:]/len_C) * (len_C+delta_cell)
    
    prim_lattice_new = np.matmul(conv_matrix_inv,lattice_new)
    
    return Structure(prim_lattice_new,structure.atomic_numbers,structure.frac_coords)

def PrepareGridPrim( struct_pmg ):
    data_path = "/work/e05/e05/isa/05_CRYSTAL/01_MnO/14_FULL_GRID_MnO/04_AF1/00_GIT_SRC/phonons/data/"
    angle_step = 0.821316026348/6
    cell_step = (4.392560947738861*0.02)/6
    conv_matrix = np.identity(3)
    num_steps = 5
    lattice_parameters = np.zeros((13,13))
    angle_parameters = np.zeros((13,13))
    supercell_matrix = np.identity(3)*2

    mno_pmg_opt = cry_gui2pmg(Crystal_gui().read_cry_gui(data_path + 'primitive/MnO_prim_af1_test.gui'))
    mno_pmg_opt_af1 = cry_gui2pmg(Crystal_gui().read_cry_gui(data_path + 'primitive/MnO_prim_af1_test.gui'))

    atoms = mno_pmg_opt_af1.atomic_numbers

    for i,delta_c in enumerate(np.arange(-num_steps,num_steps+1)):
        delta_cell = delta_c*cell_step
        #print(i,delta_cell)
        mno_pmg_tmp = set_cell(mno_pmg_opt,delta_cell,conv_matrix)
        #print(mno_pmg_tmp.lattice.a)

        for j,delta_a in enumerate(np.arange(-num_steps,num_steps+1)):
            delta_angle = delta_a*angle_step
            mno_pmg = set_angle(mno_pmg_tmp,delta_angle,conv_matrix)
            #mno_pmg_af1 = set_angle(mno_pmg_tmp,delta_angle,conv_matrix)
            #mno_pmg_af1.replace(1,1)
            #mno_pmg_af1.replace(2,1)
            mno_gui = cry_pmg2gui(mno_pmg,symmetry=False)
            a = np.round(mno_pmg.lattice.a,8)
            alpha =np.round(mno_pmg.lattice.alpha,6)
            print(delta_c, delta_a)
            print(a,alpha)
            
            mno_gui.write_crystal_gui('./MnO_cpks_%s_%s.gui'%(delta_c,delta_a),symm=False)
            mno_gui.write_crystal_gui('./MnO_freq_%s_%s.gui'%(delta_c,delta_a),symm=False)
            
            sh.copy(data_path + 'primitive/MnO_prim_af1_cpks.d12','./MnO_cpks_%s_%s.d12'%(delta_c,delta_a))
            sh.copy(data_path + 'primitive/MnO_prim_af1_freq.d12','./MnO_freq_%s_%s.d12'%(delta_c,delta_a))
            ##sh.copy('./data/primitive/MnO_prim_fm1_cpks.d12','./data/spiral/primitive/fm1/MnO_cpks_%s_%s.d12'%(delta_c,delta_a))
            #sh.copy('./data/primitive/MnO_prim_fm1_freq.d12','./data/spiral/primitive/fm1/MnO_freq_%s_%s.d12'%(delta_c,delta_a))

def PrepareGridSC( struct_pmg ):
    data_path = "/work/e05/e05/isa/05_CRYSTAL/01_MnO/14_FULL_GRID_MnO/04_AF1/00_GIT_SRC/phonons/data/"
    # SUPERCELL AF1
    angle_step = 0.821316026348/6
    cell_step = (4.392560947738861*0.02)/6
    conv_matrix = np.identity(3)
    num_steps = 5

    lattice_parameters = np.zeros((13,13))
    angle_parameters = np.zeros((13,13))

    #supercell_matrix = np.identity(3)*4
    supercell_matrix = np.identity(3)*2
    mno_pmg_opt = cry_gui2pmg(Crystal_gui().read_cry_gui(data_path + 'primitive/MnO_prim_af1_test.gui'))

    for i,delta_c in enumerate(np.arange(-num_steps,num_steps+1)):
        
        delta_cell = delta_c*cell_step
        #print(i, delta_cell)
        mno_pmg_tmp = set_cell(mno_pmg_opt,delta_cell,conv_matrix)
        
        #print(mno_pmg_tmp.lattice.a)
        for j,delta_a in enumerate(np.arange(-num_steps,num_steps+1)):

            delta_angle = delta_a*angle_step
            
            mno_pmg = set_angle(mno_pmg_tmp,delta_angle,conv_matrix)
            mno_gui = cry_pmg2gui(mno_pmg,symmetry=False)

            a = np.round(mno_pmg.lattice.a,8)
            alpha =np.round(mno_pmg.lattice.alpha,6)
            
            #mno_gui.write_crystal_gui('./data/spiral/fm1/MnO_phonons_%s_%s.gui'%(delta_c,delta_a))
            #mno_gui.write_crystal_gui('./data/spiral/fm1/MnO_phonons_%s_%s.gui'%(delta_c,delta_a))

            mno_gui.write_crystal_gui('./MnO_phonons_%s_%s.gui'%(delta_c,delta_a),symm=False)
            mno_gui.write_crystal_gui('./MnO_phonons_%s_%s.gui'%(delta_c,delta_a),symm=False)
            
            #sh.copy('./data/primitive/MnO_prim_af1_cpks.d12','./data/spiral/primitive/af1/MnO_cpks_%s_%s.d12'%(delta_c,delta_a))
            #sh.copy('./data/primitive/MnO_prim_af1_freq.d12','./data/spiral/primitive/af1/MnO_freq_%s_%s.d12'%(delta_c,delta_a))
            
            lattice_parameters[i][j] = Structure(np.matmul(conv_matrix,mno_pmg.lattice.matrix),[],[]).lattice.a
            angle_parameters[i][j] = Structure(np.matmul(conv_matrix,mno_pmg.lattice.matrix),[],[]).lattice.alpha
            
            crystal_input = Crystal_input().from_file(data_path + 'MnO_super2_af1.d12')
            #crystal_output_prim = Crystal_output().read_cry_output('./data/spiral/primitive/fm1/MnO_cpks_%s_%s.out'%(delta_c,delta_a))
            #REMOVE!!!!!!!
            #crystal_output_prim = Crystal_output().read_cry_output('./data/spiral/primitive/fm/MnO_cpks_%s_%s.out'%(delta_c,delta_a))
            #dielectric_tensor = crystal_output_prim.get_dielectric_tensor()
            crystal_input.scf_block[8] = 'ATOMSPIN\n' # tmp fix
            crystal_input.geom_block = ['EXTERNAL\n',
                                     'SCELPHONO\n',
                                     '2 0 0\n',
                                     '0 2 0\n',
                                     '0 0 2\n',
                                     'FREQCALC\n',
                                     'DISPERSI\n',
                                     'INTERPHESS\n',
                                     '10 10 10\n',
                                     '0\n',
                                     'WANG\n',
                                     '%s 0.00000  0.00000\n'%dielectric_tensor[0],
                                     '0.00000 %s 0.00000\n'%dielectric_tensor[1],
                                     '0.00000 0.00000 %s\n'%dielectric_tensor[2],
                                     'ENDFREQ\n',
                                     'END\n']
           # 
           # folder = './data/spiral/'
           # sh.copy('./data/spiral/primitive/fm/MnO_freq_%s_%s.BORN'%(delta_c,delta_a),'./data/spiral/af1/MnO_phonons_%s_%s.BORN'%(delta_c,delta_a))
            crystal_input.write_crystal_input('./MnO_phonons_%s_%s.d12'%(delta_c,delta_a))
                
def main():

    data_path = "/work/e05/e05/isa/05_CRYSTAL/01_MnO/14_FULL_GRID_MnO/04_AF1/00_GIT_SRC/phonons/data/"

    # Initial structure
    # From gui 
    structure_afm_prim_cry = Crystal_gui().read_cry_gui( ('%sMnO_afm_prim_optgeom.gui'%(data_path)) )
    #structure_prim_pmg = cry_gui2pmg(structure_prim_cry)
    mno_pmg_opt = cry_gui2pmg(structure_afm_prim_cry)

    # From pymatgen
    # Read CRYSTAL input file
    prim_crystal_input = Crystal_input().from_file(data_path + 'MnO.d12') #No GUESSP...
    scel_crystal_input = Crystal_input().from_file(data_path + 'MnO_supercell.d12') #No GUESSP...

    #Other settings
    range_lattice = [-0.01,0.01,0.005]
    range_angle = [-0.01,0.01,0.005]
    conv_matrix = np.array([[1.5, -.5, -.5],[-.5, 1.5, -.5],[-.5,-.5,1.5]]) # If working with the primitive

    # slurm file
    slurm_file = data_path + 'slurm_file.slurm'

    file = open(slurm_file, 'r')
    slurm_data = file.readlines()
    file.close()

    # Archer2
    runCRYSTAL_path = '/work/e05/e05/isa/runCRYSTAL'

    #Folders
    prim_spiral_dir = data_path + 'spiral/primitive/'
    super_spiral_dir = data_path + 'spiral/supercell/'
    material = 'MnO'

    # Correct P1 UNIT CELL #
    mno_af1 = Crystal_gui().read_cry_gui(('%sprimitive/MnO_prim_af1_test.gui'%(data_path)))
    mno_af1 = cry_gui2pmg(mno_af1)

    #PrepareGridPrim( mno_af1 )
    PrepareGridSC( mno_af1 )

if __name__ == "__main__":
   main()
