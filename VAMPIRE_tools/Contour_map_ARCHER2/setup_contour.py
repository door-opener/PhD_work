import os
import sys
import numpy as np
from shutil import copy

def getMFA( J1, J2 ):
    Tn = 2.0*(8.75)*J2
    CW = (4.0*J1)+(2.0*J2)*(8.75)
    return Tn, CW

def writeMatFile(J1, J2, fname):
    toK = np.float64( 8.61733326145E-05 )
    KtoJ = np.float64( 1.602176565E-19 )
    J1_Joules, J2_Joules = (J1*toK)*KtoJ, (J2*toK)*KtoJ
    #print("J1 %6.5fK = %s J | J2 %6.5fK = %s J" %( J1, "{:16.14E}".format( J1_Joules ), J2, "{:16.14E}".format( J2_Joules )  ) )
    loop = [ J1_Joules, J2_Joules ]
    f=open(fname, "w")
    f.write("#J1 = %6.4f J2 = %6.4f\n" %(J1,J2))
    f.write("material:num-materials=%d\n" %(3) )
    for i in range(len(loop)):
        f.write("material[%d]:material-name=Mn\n" %(i+1) )
        f.write("material[%d]:damping-constant=%3.2f\n" %(i+1, 1.0) )
        f.write("material[%d]:atomic-spin-moment=%6.4f !muB\n"%(i+1, 2.5))
        f.write('material[%d]:material-element="Mn"\n' %(i+1))
        f.write("material[%d]:density=1.00\n" %(i+1) )
        f.write("material[%d]:initial-spin-direction=random\n" %(i+1) )
        if i:
           f.write("material[%d]:exchange-matrix[%d] = %s\n" %(i+1, i+1, "{:16.14E}".format( J1_Joules ) ) )
           f.write("material[%d]:exchange-matrix[%d] = -%s\n" %(i+1, i, "{:16.14E}".format( J2_Joules ) ) )
        else:
           f.write("material[%d]:exchange-matrix[%d] = %s\n" %(i+1, i+1, "{:16.14E}".format( J1_Joules ) ) )
           f.write("material[%d]:exchange-matrix[%d] = -%s\n" %(i+1, i+2, "{:16.14E}".format( J2_Joules ) ) )
    f.write('\nmaterial[%d]:material-name="%s"\n' %(3, "O") )
    f.write('material[%d]:non-magnetic\n'%(3) )
    f.close()

def writeMCInputSTATIC(matfile, size, Tmax):
    Tinc = 2.0
    Tmin = 0.0
    numpts = Tmax/Tinc
    
    eqm_steps, loop_steps = 5000, 20000
    #print("\nTmax/min = %d/%dK | EQM steps = %d | LOOP steps = %d | #pts = %d\n" %(Tmax, Tmin, eqm_steps, loop_steps, numpts) )
    output_args = ["output:temperature", "output:material-magnetisation", "output:mean-susceptibility" ]
    pbc = [ "create:periodic-boundaries-x", "create:periodic-boundaries-y", "create:periodic-boundaries-z" ]

    cool_fn = "linear"
    precond_steps = 1000
    increment = 5

    f=open("input", "w")
    f.write("dimensions:system-size-x=%d !nm\ndimensions:system-size-y=%d !nm\ndimensions:system-size-z=%d !nm\n\n" %(size[0], size[1], size[2]))
    [ f.write("%s\n"%(x)) for x in pbc ]
    f.write("\nmaterial:file=%s\n" %(matfile) )
    f.write('material:unit-cell-file="MnO.ucf"\n')
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
    f.close()

def prepareNodeDist():

    Jst, Jend, npts = 3.0, 35.0, 32
    spacing = (Jend-Jst)/npts

    J1_range = np.zeros( npts, dtype=np.float64 )
    J2_range = np.zeros( npts, dtype=np.float64 )

    # Setting up range of J1 and J2 #
    for i in range(len(J1_range)):

        if i:
         J_cur = J_cur + spacing
        else:
         J_cur = Jst

        J1_range[i] = J_cur
        J2_range[i] = J_cur

    # Setup 2D grid #
    xx, yy = np.meshgrid( J1_range, J2_range )
    calcs_total = []

    # Set up J1 and J2 as 2D coordinates, indexed in 1D array #
    for i in range(len(xx)):
        for j in range(len(xx[i])):
            idx = ( i + (len(xx)*j ) ) + 1

            # Preparing list with element 0 as 1D index and elements 1 and 2 are J1 and J2
            current = [ idx, xx[i][j], yy[i][j] ]
            calcs_total.append( current )

    # Distribute all jobs over 8 nodes #
    node_dist = []
    for i in range( 8 ):
        cur = []
        #print("NODEID = %d" %(i) )
        for j in range(128):
            idx_1D = i + (8*j)
            cur.append( calcs_total[idx_1D] )
            #cur.append( idx_1D ) 
        node_dist.append( cur )

    return node_dist

def main():

    path_ucf = "/home/e05/e05/isa/06_VAMPIRE/00_UTILS/01_UCF_FILES/MnO.ucf"
    node_dist = prepareNodeDist() 

    size = [10, 10, 10]

    for i in range(len(node_dist)):
        print("\nNODEID = %d\n"%(i) )
        os.mkdir( "NODE%d/" %(i) )
        os.chdir( "NODE%d/" %(i) )
        for j in range(len(node_dist[i])):
            J1_cur, J2_cur = node_dist[i][j][1], node_dist[i][j][2]
            idx_1D = i + (8*j)

            os.mkdir( "STEP_%d" %( idx_1D ) )
            os.chdir( "STEP_%d" %( idx_1D ) )
            
            fname_cur = ("STEP_%d.mat" %(idx_1D))

            Tmax_cur = getMFA( J1_cur, J2_cur )[0]
            writeMatFile( J1_cur, J2_cur, fname_cur )
            writeMCInputSTATIC(fname_cur, size, Tmax_cur)

            dest = os.getcwd() + "/" + path_ucf.split("FILES/")[1]

            copy(path_ucf, dest)

            os.chdir("../")

        os.chdir("../")

if __name__=="__main__":
   main()
