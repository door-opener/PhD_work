import sys, os
import copy
import argparse
import numpy as np
from ase import Atom, Atoms, io

def writexyz(outfile, cell, pos, spec):
       spec_1, spec_2, spec_3 = list(set(spec))[0], list(set(spec))[1], list(set(spec))[2]
       s1cnt = 0
       s2cnt = 0
       s3cnt = 0
       for i in range(len(spec)):
           if spec[i] == spec_1:
              s1cnt +=1
           elif spec[i] == spec_2:
              s2cnt +=1
           else:
              s3cnt +=1
       chem_form = ("%s%d%s%d%s%d" %(spec_1, s1cnt, spec_2, s2cnt, spec_3, s3cnt))
       f=open(outfile, "w")
       f.write("%d\n" %(len(pos)) )
       f.write("%s\n" %(chem_form))
       for i in range(len(pos)):
           f.write("   %s     %4.3f     %4.3f     %4.3f\n"%(spec[i], pos[i][0], pos[i][1], pos[i][2]) )
       f.close()

def generate_tlist():
    t_list = []
    for i in range(3):
        if i > 1:
           var_i = -1
        else:
           var_i = i
        for j in range(3):
            if j > 1:
               var_j = -1
            else:
               var_j = j
            for k in range(3):
                if k > 1:
                   var_k = -1
                else:
                   var_k = k
                if i == j == k == 0:
                   continue
                else:
                   t_list.append([var_i, var_j, var_k])
    return t_list

def translateUC(t_list, ASE_cell, ASE_atoms):
 
    ASE_cell_new_atoms = Atoms(cell=[0.0,0.0,0.0])

    for i in range(len(t_list)):
        #print("Translation #%d (%d,%d,%d)" %(i, t_list[i][0], t_list[i][1], t_list[i][2]))

        for j in range(len(ASE_atoms.positions)):
            t_atom_symbol = ASE_atoms.symbols[j]
            t_atom_x = (ASE_atoms.positions[j][0] + (t_list[i][0]*ASE_cell[0][0]) )
            t_atom_y = (ASE_atoms.positions[j][1] + (t_list[i][1]*ASE_cell[1][1]) )
            t_atom_z = (ASE_atoms.positions[j][2] + (t_list[i][2]*ASE_cell[2][2]) )

            new_atom = Atom(t_atom_symbol, position=([t_atom_x, t_atom_y, t_atom_z]) )

            ASE_cell_new_atoms.append(new_atom)

    return ASE_cell_new_atoms

def translate(t_matrix, ASE_cell, ASE_pos, frac):
     if frac:
        return [ (ASE_pos[0]/ASE_cell[0][0]) + t_matrix[0], (ASE_pos[1]/ASE_cell[1][1]) + t_matrix[1], (ASE_pos[2]/ASE_cell[2][2]) + t_matrix[2] ]
     else:
        return [ ASE_pos[0] + (t_matrix[0]*ASE_cell[0][0]), ASE_pos[1] + (t_matrix[1]*ASE_cell[1][1]), ASE_pos[2] + (t_matrix[2]*ASE_cell[2][2]) ]

def distsq(ASE_a, ASE_b, ASE_cell, frac):
    distsq = np.float64(0.0)
    if frac:
       return (np.power( abs(ASE_b[0] - ASE_a[0])/ASE_cell[0][0], 2) + np.power( abs(ASE_b[1] - ASE_a[1])/ASE_cell[1][1], 2) + np.power( abs(ASE_b[2] - ASE_a[2])/ASE_cell[2][2], 2))
    else:
       return (np.power( abs(ASE_b[0] - ASE_a[0]), 2) + np.power( abs(ASE_b[1] - ASE_a[1]), 2) + np.power( abs(ASE_b[2] - ASE_a[2]), 2))

def writeExchange(num_ex_local, num_ex_ext, exchange_list_UC0, exchange_list_EXT, f, ex_type):

    #0th UC exchanges
    for i in range(num_ex_local):

        if i <= 9:
           if ex_type == "normalised-isotropic":
              f.write(" %d     %d    %d    %d    %d    %d    %10.8f\n" %(exchange_list_UC0[i][0], exchange_list_UC0[i][1], exchange_list_UC0[i][2], exchange_list_UC0[i][3], exchange_list_UC0[i][4], exchange_list_UC0[i][5], exchange_list_UC0[i][6]))
           elif ex_type == "normalised-vectorial":
              f.write(" %d     %d    %d    %d    %d    %d    %10.8f %10.8f %10.8f\n" %(exchange_list_UC0[i][0], exchange_list_UC0[i][1], exchange_list_UC0[i][2], exchange_list_UC0[i][3], exchange_list_UC0[i][4], exchange_list_UC0[i][5], exchange_list_UC0[i][6][0], exchange_list_UC0[i][6][1], exchange_list_UC0[i][6][2]))
        else:
           if ex_type == "normalised-isotropic":
              f.write("%.2d     %d    %d    %d    %d    %d    %10.8f\n" %(exchange_list_UC0[i][0], exchange_list_UC0[i][1], exchange_list_UC0[i][2], exchange_list_UC0[i][3], exchange_list_UC0[i][4], exchange_list_UC0[i][5], exchange_list_UC0[i][6]))
           elif ex_type == "normalised-vectorial":
              f.write("%.2d     %d    %d    %d    %d    %d    %10.8f %10.8f %10.8f\n" %(exchange_list_UC0[i][0], exchange_list_UC0[i][1], exchange_list_UC0[i][2], exchange_list_UC0[i][3], exchange_list_UC0[i][4], exchange_list_UC0[i][5], exchange_list_UC0[i][6][0], exchange_list_UC0[i][6][1], exchange_list_UC0[i][6][2]))

    #Translated UC exchanges
    for i in range(num_ex_ext):

        if i <= 9:
           if ex_type == "normalised-isotropic":
              f.write(" %d     %d    %d    %d    %d    %d    %10.8f\n" %(exchange_list_EXT[i][0], exchange_list_EXT[i][1], exchange_list_EXT[i][2], exchange_list_EXT[i][3], exchange_list_EXT[i][4], exchange_list_EXT[i][5], exchange_list_EXT[i][6]))
           elif ex_type == "normalised-vectorial":
              f.write(" %d     %d    %d    %d    %d    %d    %10.8f %10.8f %10.8f\n" %(exchange_list_EXT[i][0], exchange_list_EXT[i][1], exchange_list_EXT[i][2], exchange_list_EXT[i][3], exchange_list_EXT[i][4], exchange_list_EXT[i][5], exchange_list_EXT[i][6][0], exchange_list_UC0[i][6][1], exchange_list_UC0[i][6][2]))

        else:
           if ex_type == "normalised-isotropic":
              f.write("%.2d   %d    %d    %d    %d    %d   %10.8f\n" %(exchange_list_EXT[i][0], exchange_list_EXT[i][1], exchange_list_EXT[i][2], exchange_list_EXT[i][3], exchange_list_EXT[i][4], exchange_list_EXT[i][5], exchange_list_EXT[i][6]))
           elif ex_type == "normalised-vectorial":
              f.write("%.2d     %d    %d    %d    %d    %d    %10.8f %10.8f %10.8f\n" %(exchange_list_EXT[i][0], exchange_list_EXT[i][1], exchange_list_EXT[i][2], exchange_list_EXT[i][3], exchange_list_EXT[i][4], exchange_list_EXT[i][5], exchange_list_EXT[i][6][0], exchange_list_UC0[i][6][1], exchange_list_UC0[i][6][2]))

def getExch( ASE_UC0_atom, file_out, incl_J3, incl_J4, norms, cutoffs, succ ):
    ex_type = "normalised-isotropic"
    num_mat = 3

    cut_J1, cut_J2 = cutoffs[0], cutoffs[1]

    if incl_J3:
       cut_J3 = cutoffs[2]
       J3_FM_factor, J3_AFM_factor = norms[4], norms[5]
    else:
       cut_J3 = 100.0
       J3_FM_factor, J3_AFM_factor = 1.0,1.0

    if incl_J4:
       cut_J4 = cutoffs[3]
       J4_FM_factor, J4_AFM_factor = norms[6], norms[7]
    else:
       cut_J4 = 100.0
       J4_FM_factor, J4_AFM_factor = 1.0,1.0
 
    J1_FM_factor, J1_AFM_factor = norms[0], norms[1]
    J2_FM_factor, J2_AFM_factor = norms[2], norms[3]

    num_ex_local = 0
    num_ex_ext = 0

    J1_cnt_loc, J1_cnt_rem = 0, 0
    J2_cnt_loc, J2_cnt_rem = 0, 0
    J3_cnt_loc, J3_cnt_rem = 0, 0
    J4_cnt_loc, J4_cnt_rem = 0, 0

    cnt_J1_FM_loc, cnt_J1_AFM_loc = 0, 0
    cnt_J2_FM_loc, cnt_J2_AFM_loc = 0, 0
    cnt_J3_FM_loc, cnt_J3_AFM_loc = 0, 0
    cnt_J4_FM_loc, cnt_J4_AFM_loc = 0, 0

    cnt_J1_FM_rem, cnt_J1_AFM_rem = 0, 0
    cnt_J2_FM_rem, cnt_J2_AFM_rem = 0, 0
    cnt_J3_FM_rem, cnt_J3_AFM_rem = 0, 0
    cnt_J4_FM_rem, cnt_J4_AFM_rem = 0, 0

    #Generate translation list.
    cell = ASE_UC0_atom.get_cell()
    t_list = generate_tlist()
    uc_size = [ cell[0][0], cell[1][1], cell[2][2] ]

    exchange_list_UC0 = []
    exchange_list_EXT = []

    excl = "O"

    #Exchange list 0th UC.
    for i in range(len(ASE_UC0_atom.positions)):

        if ASE_UC0_atom.symbols[i] != excl:
           #NOTE PURPOSEFULLY DID NOT USE N-1! BECAUSE INTERACTIONS ARE NOT SYMMETRIC!
           for j in range(len(ASE_UC0_atom.positions)):

               if ASE_UC0_atom.symbols[j] != excl:
                  if i != j:
                     dist = distsq( ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j], cell, frac=False)

                     if cut_J1 < dist < cut_J2:
                        if ASE_UC0_atom.symbols[i] != ASE_UC0_atom.symbols[j]:
                            cnt_J2_AFM_loc += 1
                            exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J2_AFM_factor] )
                        else:
                            cnt_J2_FM_loc += 1
                            exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J2_FM_factor] )

                        num_ex_local += 1
                        J2_cnt_loc += 1

                     elif dist < cut_J1:
                        if ASE_UC0_atom.symbols[i] != ASE_UC0_atom.symbols[j]:
                            cnt_J1_AFM_loc += 1
                            #print("i (%d) typ = %s | j typ = %s" %( i, ASE_UC0_atom.symbols[i], ASE_UC0_atom.symbols[j] ) )
                            #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                            #print("J1 (%d, %d, AFM) dist = %4.3f\n" %(i, j, np.sqrt(dist)) )
                            exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J1_AFM_factor] )
                        else:
                            cnt_J1_FM_loc += 1
                            #print("i (%d) typ = %s | j typ = %s" %( i, ASE_UC0_atom.symbols[i], ASE_UC0_atom.symbols[j] ) )
                            #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                            #print("J1 (%d, %d, FM) dist = %4.3f\n" %(i, j, np.sqrt(dist)) )
                            exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J1_FM_factor] )

                        num_ex_local += 1
                        J1_cnt_loc += 1

                     elif cut_J2 < dist < cut_J3:
                        if incl_J3:
                           if ASE_UC0_atom.symbols[i] != ASE_UC0_atom.symbols[j]:
                               cnt_J3_AFM_loc += 1
                               exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J3_AFM_factor] )
                           else:
                               cnt_J3_FM_loc += 1
                               exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J3_FM_factor] )
                          #print("i (%d) typ = %s | j typ = %s" %( i, ASE_UC0_atom.symbols[i], ASE_UC0_atom.symbols[j] ) )
                          #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                          #print("dist(i, j) = %4.3f\n" %(np.sqrt(dist)) )
                           num_ex_local += 1
                           J3_cnt_loc += 1

                     elif cut_J3 < dist < cut_J4:
                        if incl_J4:
                           if ASE_UC0_atom.symbols[i] != ASE_UC0_atom.symbols[j]:
                               cnt_J4_AFM_loc += 1
                               exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J4_AFM_factor] )
                           else:
                               cnt_J4_FM_loc += 1
                               exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0, J4_FM_factor] )
                          #print("i (%d) typ = %s | j typ = %s" %( i, ASE_UC0_atom.symbols[i], ASE_UC0_atom.symbols[j] ) )
                          #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                          #print("dist(i, j) = %4.3f\n" %(np.sqrt(dist)) )
                           num_ex_local += 1
                           J4_cnt_loc += 1

    #Prepare Global list from translations of UC0.
    #If Using cartesian cutoffs, multiply the translation by the cell vectors.
    for i in range(len(t_list)):
        local = []
        remote = []

        #Translate the whole unit cell first
        for j in range(len(ASE_UC0_atom.positions)):
            untrans = [ ASE_UC0_atom.positions[j][0], ASE_UC0_atom.positions[j][1], ASE_UC0_atom.positions[j][2] ]
            local.append( untrans )

            trans = translate( t_list[i], cell, ASE_UC0_atom.positions[j], frac=False )
            remote.append( trans )
           #print("Translation #%d [ %d %d %d ]\n" %(i, t_list[i][0], t_list[i][1], t_list[i][2]))

        for k in range( len(local) ):
            if ASE_UC0_atom.symbols[k] != excl:
               for l in range( len(remote) ):
                   if ASE_UC0_atom.symbols[l] != excl:
                      dist = distsq( local[k], remote[l], cell, frac=False )

                      if cut_J1 < dist < cut_J2:
                         if ASE_UC0_atom.symbols[k] != ASE_UC0_atom.symbols[l]:
                             cnt_J2_AFM_rem += 1
                             exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J2_AFM_factor] )
                         else:
                             cnt_J2_FM_rem += 1
                             exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J2_FM_factor] )
                             #print("AFM\n")
                             #print("i (%d) typ = %s | j typ = %s" %( i, ASE_UC0_atom.symbols[k], ASE_UC0_atom.symbols[l] ) )
                             #print("J3 Local (i=%d) %s, Remote (j=%d) %s" %(k, local[k], l, remote[l]))
                             #print("DIST = %4.3f\n" %(np.sqrt(dist) ))

                         num_ex_ext += 1
                         J2_cnt_rem += 1

                      elif dist < cut_J1:
                         if ASE_UC0_atom.symbols[k] != ASE_UC0_atom.symbols[l]:
                             exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J1_AFM_factor] )
                             cnt_J1_AFM_rem += 1
                         else:
                             exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J1_FM_factor] )
                             cnt_J1_FM_rem += 1
                         #print("J2 Local (i=%d) %s, Remote (j=%d) %s" %(k, local[k], l, remote[l]))
                         #print("DIST = %4.3f\n" %(dist) )
                         #print("J2 (%d, %d) dist = %4.3f\n" %(i, j, dist) )
                         #exchange_list_UC0.append( [num_ex_local, i, j, 0, 0, 0 ] )
                         num_ex_ext += 1
                         J1_cnt_rem += 1

                      elif cut_J2 < dist < cut_J3:
                         if incl_J3:
                            if ASE_UC0_atom.symbols[k] != ASE_UC0_atom.symbols[l]:
                                cnt_J3_AFM_rem += 1
                                exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J3_AFM_factor ] )
                            else:
                                cnt_J3_FM_rem += 1
                                exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J3_FM_factor ] )
                           #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                           #print("dist(i, j) = %4.3f\n" %(np.sqrt(dist)) )
                            num_ex_ext += 1
                            J3_cnt_rem += 1

                      elif cut_J3 < dist < cut_J4:
                         if incl_J4:
                            if ASE_UC0_atom.symbols[k] != ASE_UC0_atom.symbols[l]:
                                cnt_J4_AFM_rem += 1
                                exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J4_AFM_factor ] )
                            else:
                                cnt_J4_FM_rem += 1
                                exchange_list_EXT.append( [ num_ex_local + num_ex_ext, k, l, t_list[i][0], t_list[i][1], t_list[i][2], J4_FM_factor ] )
                          #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                          #print("i = %s j = %s"%(ASE_UC0_atom.positions[i], ASE_UC0_atom.positions[j]))
                          #print("dist(i, j) = %4.3f\n" %(np.sqrt(dist)) )
                            num_ex_ext += 1
                            J4_cnt_rem += 1

    if succ == "info":
       print("         a = %3.2f, b = %3.2f, c = %3.2f" %(uc_size[0], uc_size[1], uc_size[2]) )
       print("---------------------------------------------------")
       print("# Sublattices                    = %d" %( num_mat ) )
       print("---------------------------------------------------")
       print("Exchange Info:")
       print("---------------------------------------------------")
       print("Exch Type: normalised-isotropic\n")
       print(" J1 cutoff                         = %4.3f" %(np.sqrt(cut_J1)) )
       print(" J2 cutoff                         = %4.3f" %(np.sqrt(cut_J2)) )
       print(" J3 cutoff                         = %4.3f" %(np.sqrt(cut_J3)) )
       print(" J4 cutoff                         = %4.3f" %(np.sqrt(cut_J4)) )
       print("---------------------------------------------------")
       print("Total Exchanges In Unit Cell 0      = %d" %(num_ex_local) )
       print("---------------------------------------------------")
       print("  #J1                               = %d" %(J1_cnt_loc ) )
       print("    #J1 FM                          = %d" %(cnt_J1_FM_loc))
       print("    #J1 AFM                         = %d" %(cnt_J1_AFM_loc))
       print("\n  #J2                               = %d" %(J2_cnt_loc ) )
       print("    #J2 FM                          = %d" %(cnt_J2_FM_loc))
       print("    #J2 AFM                         = %d" %(cnt_J2_AFM_loc))
       print("\n  #J3                               = %d" %(J3_cnt_loc ) )
       print("    #J3 FM                          = %d" %(cnt_J3_FM_loc))
       print("    #J3 AFM                         = %d" %(cnt_J3_AFM_loc))
       print("\n  #J4                               = %d" %(J4_cnt_loc ) )
       print("    #J4 FM                          = %d" %(cnt_J4_FM_loc))
       print("    #J4 AFM                         = %d" %(cnt_J4_AFM_loc))
       print("---------------------------------------------------")
       print("Total Exchanges In Translated Cells = %d" %(num_ex_ext) )
       print("---------------------------------------------------")
       print("  #J1                               = %d" %(J1_cnt_rem ) )
       print("    #J1 FM                          = %d" %(cnt_J1_FM_rem))
       print("    #J1 AFM                         = %d" %(cnt_J1_AFM_rem))
       print("\n  #J2                               = %d" %(J2_cnt_rem ) )
       print("    #J2 FM                          = %d" %(cnt_J2_FM_rem))
       print("    #J2 AFM                         = %d" %(cnt_J2_AFM_rem))
       print("\n  #J3                               = %d" %(J3_cnt_rem ) )
       print("    #J3 FM                          = %d" %(cnt_J3_FM_rem))
       print("    #J3 AFM                         = %d" %(cnt_J3_AFM_rem))
       print("\n  #J4                               = %d" %(J4_cnt_rem ) )
       print("    #J4 FM                          = %d" %(cnt_J4_FM_rem))
       print("    #J4 AFM                         = %d" %(cnt_J4_AFM_rem))
       print("---------------------------------------------------")
       print("Total Exchanges                     = %d" %(num_ex_local + num_ex_ext) )
       print("---------------------------------------------------")
       sys.exit()

    f=open(file_out, "w")
    f.write("#UCSIZE\n")
    f.write("%3.2f    %3.2f    %3.2f\n" %(uc_size[0], uc_size[1], uc_size[2]) )
    f.write("#CELLVECTORS\n")
    f.write("%3.2f    %3.2f    %3.2f\n" %(1.0, 0.0, 0.0) )
    f.write("%3.2f    %3.2f    %3.2f\n" %(0.0, 1.0, 0.0) )
    f.write("%3.2f    %3.2f    %3.2f\n" %(0.0, 0.0, 0.0) )
    f.write("#ATOMINFO\n")
    f.write("%d    %d\n" %( len(ASE_UC0_atom), num_mat) )

    #UC_info
    for i in range(len(ASE_UC0_atom.positions)):
        mat_id = None

        if ASE_UC0_atom.symbols[i] == "Ca":
           mat_id = 0
        elif ASE_UC0_atom.symbols[i] == "Ba":
           mat_id = 1
        else:
           mat_id = 2

        if i <= 9:
           f.write("%d    %4.3f    %4.3f    %4.3f    %d    %d    %d\n" %(i, ASE_UC0_atom.positions[i][0]/cell[0][0], ASE_UC0_atom.positions[i][1]/cell[1][1], ASE_UC0_atom.positions[i][2]/cell[2][2], mat_id, 0, 0) )
        else:
           f.write("%.2d    %4.3f    %4.3f    %4.3f    %d    %d    %d\n" %(i, ASE_UC0_atom.positions[i][0]/cell[0][0], ASE_UC0_atom.positions[i][1]/cell[1][1], ASE_UC0_atom.positions[i][2]/cell[2][2], mat_id, 0, 0) )

    f.write("#INTERACTION_LIST\n")
    f.write("%d    %s\n" %(num_ex_local + num_ex_ext, ex_type) )
    writeExchange(num_ex_local, num_ex_ext, exchange_list_UC0, exchange_list_EXT, f, ex_type)
    f.close()
    print("File written as [%s]" %(file_out) )

def readInput():
    f=open("UCF.in", "r")
    for line in f:
        split = line.split("=")
        if "geometry" in line:
            geom_file=(split[-1]).strip()
        if "NN_range" in line:
            NN_range=int(split[-1])
        if "cutoffs" in line:
            cutoffs = (split[-1].strip()).split(",")
        if "norms" in line:
            norms = (split[-1].strip()).split(",")
    f.close()
    return geom_file, NN_range, np.asarray(cutoffs, dtype=np.float64), np.asarray(norms, dtype=np.float64)    

def main():
 
    """
    This code is used to generate .ucf files for use with the magnetic MC software VAMPIRE.
    
    Atomic structure must be provided as a .cif file, the filename is given in the input file 'UCF.in' under geometry='<input-file.cif>'.

    The number of nearest neighbours is given in 'UCF.in' under NN_range. For first and second nearest neighbours only, set NN_range to 2.

    Distance cutoffs are also supplied in UCF.in, under cutoffs. The number of values given for cutoffs should be equal to the NN_range.
    NB. Distance cutoffs must be squared when provided in UCF.in. i.e. if your cutoff is 2 Angs, set the cutoff to 4 Angs.

    Normalisation constants are included in UCF.in, under norms.

    Normalisation constants can be used to adjust the strength of magnetic coupling for specific interactions.

      This is necessary if the number of exchange interactions you want to include exceeds the number of materials specified in the material file (.mat).
 
      These can also be used to adjust the strength of interactions between spins within the same cutoff range.
      For example, when you have pairs of spins in a given NN shell who's interaction strengths differ.

    """

    # To test translations #
    #t_list = generate_tlist()
    #TEST = translateUC(t_list, cell, ASE_UC0_atom)
    #spec = TEST.symbols
    #pos = TEST.positions
    #cell = [ cell[0][0], cell[1][1], cell[2][2] ]
    
    if os.path.exists("UCF.in") == False:    
       print("No input file found.")
       print("Exiting...")
       sys.exit()

    else:
       geom_file, NN_range, cutoffs, norms = readInput()
       ASE_UC0_atom = io.read(geom_file)

       print("---------------------------------------------------")
       print("                  UCF BUILDER                      ")
       print("---------------------------------------------------")
       print("Input geometry = %s" %(geom_file) )
       print("# of Nearest Neighbours = %d" %(NN_range) )
       print("---------------------------------------------------")

       try:
          option = sys.argv[1]

       except IndexError:
          print("No option entered!")
          print("Exiting...")
          sys.exit()

       incl_J3, incl_J4 = False, False

       if NN_range > 4:
          print("Maximum of four nearest neighbours.") 
          print("Exiting...")
          sys.exit()

       elif NN_range == 3:
          incl_J3 = True

       elif NN_range == 4:
          incl_J3 = True
          incl_J4 = True

       if option == "-i":
          getExch( ASE_UC0_atom, "dummy.out", incl_J3, incl_J4, norms, cutoffs, "info" )

       elif option == "-b": 
          getExch( ASE_UC0_atom, "UCF.out", incl_J3, incl_J4, norms, cutoffs, "build" )

       else:
          print("Unrecognised option entered!")
          print("Exiting...")
          sys.exit()

if __name__=="__main__":
   main()
