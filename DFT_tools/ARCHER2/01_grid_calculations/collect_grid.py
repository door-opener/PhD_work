import sys
import os
import numpy as np
import subprocess

def pyGrep(grep_args):
    """
    grep_args: [0] - grep string.
               [1] - file to be 'grepped'.
               [2] - grep arguement modifiers.
               [3] - argument for PIPE command.
    """
    #standard grep
    if len(grep_args) == 2:
        gstring = ("grep '%s' %s" %(str(grep_args[0]), str(grep_args[1]) ) )
    #if args == 3, we are using a grep argument.
    elif len(grep_args) == 3:
        gstring = ("grep %s '%s' %s" %(str(grep_args[2]), str(grep_args[0]), str(grep_args[1]) ) )
    #PIPE and grep.
    else:
        try:
            gstring = ("grep %s '%s' %s | %s" %(grep_args[2], grep_args[0], grep_args[1], grep_args[3]))
        except:
           print("too many args!\n")
    out = subprocess.check_output(gstring, stderr=subprocess.STDOUT, shell=True).decode('utf-8')
    return str(out)

def getSCF(filename):
    ener, steps = pyGrep( [ "== SCF ENDED", filename ] ).split()[8], pyGrep( [ "== SCF ENDED", filename ] ).split()[-1]
    return np.float64(ener), int(steps)

def getOuts( in_dir ):
    shit = sorted( [ x for x in os.listdir( in_dir ) if ".sh" not in x ], key= lambda x: int(x.split("_")[0]) )
    return [ ("%s/%s/MnO.out"%(in_dir,x)) for x in shit ]

def getData( out_cur ):
    toeV = 27.211384500
    try:
       eg_up, eg_down = np.float64(pyGrep( ["ALPHA      ELECTRONS", out_cur, "-A5", "tail -1" ] ).split()[4]), np.float64(pyGrep( ["BETA       ELECTRONS", out_cur, "-A5", "tail -1" ] ).split()[4])
    except:
       eg_up, eg_down = np.float64(pyGrep( ["ALPHA      ELECTRONS", out_cur, "-A4", "tail -1" ] ).split()[4]), np.float64(pyGrep( ["BETA       ELECTRONS", out_cur, "-A4", "tail -1" ] ).split()[4])
    sd = np.float64( pyGrep( ["TOTAL ATOMIC SPINS  ", out_cur, "-A1", "tail -1"] ).split()[0] )
    ener, steps = getSCF( out_cur )
    time = np.float64( pyGrep([ "END         TELAPSE", out_cur ]).split()[3] )
    return ener*toeV, steps, sd, [ eg_up, eg_down ], time

def main():

    kb = 0.00008617333262

    verbose = True

    try:
      in_dir = sys.argv[1]
    except:
      print("\nEnter a directory.\n")
      sys.exit()

    AFII_src = [ ("%s%s" %(in_dir, x)) for x in os.listdir( in_dir ) if "AFII" in x ][0]
    #AFI_src = [ ("%s%s" %(in_dir, x)) for x in os.listdir( in_dir ) if "AFI" in x and "AFII" not in x ][0]
    FM_src = [ ("%s%s" %(in_dir, x)) for x in os.listdir( in_dir ) if "FM" in x ][0]

    AFII_outs = getOuts( AFII_src )
    FM_outs = getOuts( FM_src )

    AFII_write, FM_write = [], []
    EFMAFII_write = []

    for i in range(len(AFII_outs)):

        #lc_cur = np.float64(( AFII_outs[i].split("_")[2] ).replace("A",""))
        lc_cur = np.float64(( AFII_outs[i].split("_")[3] ).replace("A",""))
        #lc_cur = np.float64(( AFII_outs[i].split("_")[5] ).replace("A",""))
        ang_cur = np.float64( ( AFII_outs[i].split("_")[-1]).replace("ALP/MnO.out",""))

        AFII_E, AFII_st, AFII_sd, AFII_eg, AFII_time = getData( AFII_outs[i] )
        FM_E, FM_st, FM_sd, FM_eg, FM_time = getData( FM_outs[i] )

        AFII_write.append( [lc_cur, ang_cur, AFII_E, AFII_eg, AFII_sd/2.0, AFII_st, AFII_time ] )
        #AFI_write.append( [lc_cur, ang_cur, AFI_E/8.0, AFI_eg, AFI_sd/4.0, AFI_st, AFI_time ] )
        FM_write.append( [lc_cur, ang_cur, FM_E, FM_eg, FM_sd/2.0, FM_st, FM_time ] )
        EFM_AFII_cur = (AFII_E - FM_E)
        EFMAFII_write.append( [lc_cur, ang_cur, EFM_AFII_cur] )

    if verbose:

       print("\n-------------------------------------------------------------------\n")
       print(" a0     alpha       EAFII (eV)      Eg alpha (eV)    Eg beta (eV)    S/2 (ub)  #SCF   Time(s) ")
       print("\n-------------------------------------------------------------------\n")
       [ print("%3.2f    %3.2f    %16.14f    %6.4f    %6.4f   %6.4f   %2d    %6.4f"%(x[0], x[1], x[2], x[3][0], x[3][1], x[4], x[5], x[6])) for x in AFII_write ]

       #print("\n-------------------------------------------------------------------\n")
       #print(" a0     alpha       EAFI (eV/atom)      Eg (eV)    S/2 (ub)  #SCF   Time(s) ")
       #print("\n-------------------------------------------------------------------\n")
       #[ print("%3.2f    %3.2f    %16.14f    %6.4f    %6.4f   %2d    %6.4f"%(x[0], x[1], x[2], x[3], x[4], x[5], x[6])) for x in AFI_write ]

       print("\n-------------------------------------------------------------------\n")
       print(" a0     alpha       EFM (eV)      Eg alpha (eV)    Eg beta (eV)    S/2 (ub)  #SCF   Time(s) ")
       print("\n-------------------------------------------------------------------\n")
       [ print("%3.2f    %3.2f    %16.14f    %6.4f    %6.4f   %6.4f   %2d    %6.4f"%(x[0], x[1], x[2], x[3][0], x[3][1], x[4], x[5], x[6])) for x in FM_write ]
     

       print("\n-------------------------------------------------------------------\n")
       print(" a0    alpha    EAFII - EFM (eV)    ")
       print("\n-------------------------------------------------------------------\n")
       [ print("%3.2f    %3.2f   %16.14f" %(x[0], x[1], x[2])) for x in EFMAFII_write ]

       #print("\n-------------------------------------------------------------------\n")
       #print(" a0    alpha    EFM - EAFI (eV/atom)    ")
       #print("\n-------------------------------------------------------------------\n")
       #[ print("%3.2f    %3.2f   %16.14f" %(x[0], x[1], x[2])) for x in EFMAFI_write ]

       #print("\n-------------------------------------------------------------------\n")
       #print(" a0    alpha    J1 (eV/atom)    J1 (meV/atom)    J1 (K/atom)    ")
       #print("\n-------------------------------------------------------------------\n")
       #[ print("%3.2f    %3.2f   %12.10f    %12.10f    %6.4f" %(x[0], x[1], x[2], x[2]*1000.0, x[2]/kb)) for x in J1_write ]
 
       #print("\n-------------------------------------------------------------------\n")
       #print(" a0    alpha    J2 (eV/atom)    J2 (meV/atom)    J2 (K/atom)    ")
       #print("\n-------------------------------------------------------------------\n")
       #[ print("%3.2f    %3.2f   %12.10f    %12.10f    %6.4f" %(x[0], x[1], x[2], x[2]*1000.0, x[2]/kb)) for x in J2_write ]

if __name__ == "__main__":
   main()
