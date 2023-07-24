import sys
import os
import numpy as np
from numpy.polynomial import Polynomial

def ReadDat( in_file ):
    a0 = []
    E_FM, E_AF1, E_AF2, S = [], [], [], []
    f=open( in_file )
    next(f)
    for line in f:
        split = line.split()
        a0.append(split[0])
        E_FM.append(split[1])
        E_AF1.append(split[4])
        E_AF2.append(split[7])
        S.append(split[-1])
    f.close()
    return np.asarray(a0, dtype=np.float64), np.asarray(E_FM, dtype=np.float64), np.asarray(E_AF1, dtype=np.float64), np.asarray(E_AF2, dtype=np.float64), np.asarray(S, dtype=np.float64)

def Fit( x,y,order, mina=True ):
    fitted, lsq = Polynomial.fit(x, y, deg=order, full=True)
    min_bounds = [ min(x), max(x) ]
    #print("--------------------------")
    #[ print("x^%d = %s" %(i, "{:16.14E}".format(x) )) for i,x in enumerate(fitted.convert().coef) ]
    #print("--------------------------")
    #print("Rsq = %s" %("{:16.14E}".format(lsq[0][0]) ) )
    #print("--------------------------")
    if mina:
       eqm_a0 = [ x for x in fitted.deriv().roots() if x.imag == 0.0 and min_bounds[0] < x < min_bounds[1] ][0]
       return eqm_a0, fitted, lsq
    else:
       return fitted, lsq

def GetOptJs( FM, AF1, AF2, S, a0, a0_min, NiO=False):
    kB = 0.00008617333262
    toJoules = np.float64( 1.602176565E-19 )
    dummy1, fit_AF1, lsq_AF1 = Fit( a0, AF1, 7 )
    dummy2, fit_AF2, lsq_AF2 = Fit( a0, AF2, 7 )
    dummy3, fit_FM, lsq_FM = Fit( a0, FM, 7 )

    fit_S, lsq_S = Fit( a0, S, 7, mina=False )

    coef1 = fit_AF1.convert().coef
    coef2 = fit_AF2.convert().coef
    coef3 = fit_FM.convert().coef
    coef4 = fit_S.convert().coef
    E1, E2, E3, S_opt, s_sq = 0.0, 0.0, 0.0, 0.0, 0.0

    for i in range(len(coef1)):
        E1 += np.power( a0_min, i )*coef1[i]
        E2 += np.power( a0_min, i )*coef2[i]
        E3 += np.power( a0_min, i )*coef3[i]
        S_opt += np.power( a0_min, i )*coef4[i]

    if NiO:
        s_sq = 1.0
        s_dft = np.power(S_opt/4,2)
    else:
        s_sq = 6.25
        s_dft = np.power(S_opt/4,2)

    print("------------------------------")
    print("Optimal a0 = %4.3f" %(a0_min))
    print("------------------------------")
    print("AF1 = " + str(E1) )
    print("AF2 = " + str(E2) )
    print("FM  = " + str(E3) )
    print("S   = " + str(S_opt/4.0) + " uB")
    print("------------------------------")
    print("\nFM - AF1 = %10.8f meV/fu" %((E3-E1)*1000.0) )
    print("FM - AF2 = %10.8f meV/fu\n" %((E3-E2)*1000.0) )
    print("------------------------------")

    # Equations from Logsdail et al. #
    J1_opt = (E3-E1)/(8.0*s_sq)
    J2_opt = (E3-E2)/(6.0*s_sq) - J1_opt

    J1_opt1 = (E3-E1)/(8.0*np.power(S_opt/4,2))
    J2_opt1 = (E3-E2)/(6.0*np.power(S_opt/4,2)) - J1_opt1

    # Energy correction #
    #J2_opt = J2_opt*0.9
    #J2_opt1 = J2_opt1*0.9

    print("------------------------------")
    print("S = 5/2")
    print("------------------------------")
    print("J1 = %s eV | %8.6f meV | %6.4f K" %("{:14.12E}".format( J1_opt ), J1_opt*1000.0, J1_opt/kB  ) )
    print("J2 = %s eV | %8.6f meV | %6.4f K" %("{:14.12E}".format(J2_opt), J2_opt*1000.0, J2_opt/kB ) )
    print("------------------------------")
    print("J1 = %s J " %("{:14.12E}".format( J1_opt*toJoules ) ) )
    print("J2 = %s J " %("{:14.12E}".format(J2_opt*toJoules ) ) )
    print("------------------------------")
    print("TN (MFA)       = %8.6f K" %( 2.0*s_sq*(J2_opt/kB) ) )
    print("TN (MFA corr. [Ashcroft & Mermin] ) = %8.6f K" %( ( 2.0*s_sq*(J2_opt/kB) )*0.75  ) )
    print("TN (MFA corr. [Garnin et al.]     ) = %8.6f K" %( ( 2.0*s_sq*(J2_opt/kB) )*0.79 ) )
    print("CW (MFA)       = %8.6f K" %( s_sq*( (4.0*J1_opt/kB) + (2.0*J2_opt/kB) ) ) )
    print("------------------------------")

    print("------------------------------")
    print("S = S_dft/2")
    print("------------------------------")
    print("J1 = %s eV | %8.6f meV | %6.4f K" %("{:14.12E}".format( J1_opt1 ), J1_opt1*1000.0, J1_opt1/kB  ) )
    print("J2 = %s eV | %8.6f meV | %6.4f K" %("{:14.12E}".format(J2_opt1), J2_opt1*1000.0, J2_opt1/kB ) )
    print("------------------------------")
    print("J1 = %s J " %("{:14.12E}".format( J1_opt1*toJoules ) ) )
    print("J2 = %s J " %("{:14.12E}".format(J2_opt1*toJoules ) ) )
    print("------------------------------")
    print("TN (MFA)       = %8.6f K" %( 2.0*s_dft*(J2_opt1/kB) ) )
    print("TN (MFA corr. [Ashcroft & Mermin] ) = %8.6f K" %( ( 2.0*s_dft*(J2_opt1/kB) )*0.75  ) )
    print("TN (MFA corr. [Garnin et al.]     ) = %8.6f K" %( ( 2.0*s_dft*(J2_opt1/kB) )*0.79 ) )
    print("CW (MFA)       = %8.6f K" %( 6.25*( (4.0*J1_opt1/kB) + (2.0*J2_opt1/kB) ) ) )
    print("------------------------------\n")

    if NiO:
       return [ J1_opt*toJoules*-1.0, J2_opt*toJoules*-1.0, 1.00 ], [ J1_opt1*toJoules*-1.0, J2_opt1*toJoules*-1.0, S_opt/4.0 ]
    else:
       return [ J1_opt*toJoules*-1.0, (J2_opt)*toJoules*-1.0, 2.50 ], [ J1_opt1*toJoules*-1.0, (J2_opt1)*toJoules*-1.0, S_opt/4.0 ]

def writeMat( S_formal, S_DFT, in_file, NiO=False, hyb=False ):
    dump = [ S_formal[-1], S_DFT[-1] ]
    out_file = in_file.replace("_in.dat", ".mat")

    fmtJ = "{:12.10E}"

    for i in range(len( dump)):
        if i == 0:
           out_file = in_file.replace("_in.dat", "_Sformal.mat")
           J1, J2 = S_formal[0], S_formal[1]
        else:
           out_file = in_file.replace("_in.dat", "_SDFT.mat")
           J1, J2 = S_DFT[0], S_DFT[1]

        f=open( out_file, "w")
        f.write("##############################################\n")
        if hyb:
           f.write("#%s, J2/J1 ~ %1.2f\n" %(in_file.replace("_in.dat", ""),(J2*0.9)/J1) )
           f.write("#J2 HAS BEEN CORRECTED!\n")
        else:
           f.write("#%s, J2/J1 ~ %1.2f\n" %(in_file.replace("_in.dat", ""),(J2)/J1) )
        f.write("##############################################\n")
        f.write("#DAMPING!\n")
        f.write("#If LLG-HEUN -> 0.1 - 1.0\n")
        f.write("#If MC -> 1.0\n")
        f.write("##############################################\n")
        f.write("material:num-materials=3\n")
        f.write("##############################################\n")

        f.write("###################MAT1#######################\n")
        if NiO:
           f.write("material[1]:material-name=Ni\n")
        else:
           f.write("material[1]:material-name=Mn\n")
        f.write("material[1]:damping-constant=1.00\n")
        f.write("material[1]:atomic-spin-moment=%6.4f !muB\n" %(dump[i]))
        if NiO:
           f.write('material[1]:material-element="Ni"\n')
        else:
           f.write('material[1]:material-element="Mn"\n')
        f.write("material[1]:density=1.00\n")
        f.write("material[1]:initial-spin-direction=random\n")
        f.write("material[1]:exchange-matrix[1] =  %s\n" %( fmtJ.format(J1) ) )
        if hyb:
           f.write("material[1]:exchange-matrix[2] =  %s\n" %( fmtJ.format(J2*0.9) ) )
        else:
           f.write("material[1]:exchange-matrix[2] =  %s\n" %( fmtJ.format(J2) ) )

        f.write("###################MAT2#######################\n")
        if NiO:
           f.write("material[2]:material-name=Ni\n")
        else:
           f.write("material[2]:material-name=Mn\n")
        f.write("material[2]:damping-constant=1.00\n")
        f.write("material[2]:atomic-spin-moment=%6.4f !muB\n" %(dump[i]))
        if NiO:
           f.write('material[2]:material-element="Ni"\n')
        else:
           f.write('material[2]:material-element="Mn"\n')
        f.write("material[2]:density=1.00\n")
        f.write("material[2]:initial-spin-direction=random\n")
        f.write("material[2]:exchange-matrix[2] =  %s\n" %( fmtJ.format(J1) ) )
        if hyb:
           f.write("material[2]:exchange-matrix[1] =  %s\n" %( fmtJ.format(J2*0.9) ) )
        else:
           f.write("material[2]:exchange-matrix[1] =  %s\n" %( fmtJ.format(J2) ) )

        f.write("###################MAT3#######################\n")
        f.write('material[3]:material-name="O"\n')
        f.write("material[3]:non-magnetic\n")
        f.close()
    
def main():

   write = False

   try:
      in_dir = sys.argv[1]
      dummy = [ x for x in os.listdir( in_dir ) if ".dat" in x ]
      dirs = [ ("%s/%s" %(in_dir, x)) for x in dummy ]
   except:
     print("\nEnter a directory.\n")
     sys.exit()

   for i in range(len( dirs )):
       print("FILE CUR = %s\n" %(dirs[i]) )
       order = 7
       a0, FM, AF1, AF2, S = ReadDat( dirs[i] )

       mina0 = Fit( a0, AF2, order )
       print("min a0 AF2 = %8.6f A" %(mina0[0].real))

       if "NiO" in dirs[i]:
          S_formal, S_DFT = GetOptJs( FM, AF1, AF2, S, a0, mina0[0].real, NiO=True )
          if write:
             writeMat( S_formal, S_DFT, dirs[i], NiO=True )
       else:
          S_formal, S_DFT = GetOptJs( FM, AF1, AF2, S, a0, mina0[0].real )
          if write:
             if "PBESol0" in dirs[i] or "PBE0" in dirs[i]:
                writeMat( S_formal, S_DFT, dirs[i], NiO=False, hyb=True)
             else:
                writeMat( S_formal, S_DFT, dirs[i], NiO=False, hyb=False )

   sys.exit()

   order = 7
   a0, FM, AF1, AF2, S = ReadDat( in_file )
   mina0 = Fit( a0, AF2, order )
   print("min a0 FM = %8.6f A" %(mina0[0].real))
   GetOptJs( FM, AF1, AF2, S, a0, mina0[0].real )

if __name__ == "__main__":
   main()
