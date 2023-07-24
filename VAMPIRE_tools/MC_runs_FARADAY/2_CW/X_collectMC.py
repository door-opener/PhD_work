import os
import sys
import numpy as np
import subprocess as sp

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
    out = sp.check_output(gstring, stderr=sp.STDOUT, shell=True).decode('utf-8')
    return str(out)

def getMolX( infile, X ):
    toM3 = 1e+27
    grep_str = "dimensions:system-size-x"
    dummy = pyGrep([grep_str, infile]).split("=")[-1].replace("!nm", "")
    volnm = np.power( np.float64(dummy), 3)/toM3
    #density = (5.37/1000.0)/1e+6 # kg/m^3
    mol_mass = 54.938049 # g/mol
    density = 5.37       # g/cm^3
    #mol_mass = 54.938049/1000.0 # kg/mol
    #print(density, mol_mass)
    return (mol_mass/density)*X

def main():

    plot = False

    try:
      dirs = sorted([ ("%s%s"%(sys.argv[1], x)) for x in os.listdir( sys.argv[1] ) ])
    except:
      print("Please enter directory!")
      sys.exit()

    files = [ ("%s/%s/output" %(os.getcwd(), x)) for x in dirs ]

    print("\n------------------------------------------------------------------------------\n")

    print("T    X    X^-1")

    fit_x = np.zeros(len(files), dtype=np.float64 )
    fit_y = np.zeros(len(files), dtype=np.float64 )

    for i in range(len(files)):

        print(files[i])
        continue

        X_xyz, X = [], []
        T = np.float64(0.0)

        f=open(files[i], "r")
        [ next(f) for x in range(8) ]
        for j,line in enumerate(f):
            split = line.split()

            if j == 0:
               T = np.float64(split[0] )

            X_xyz.append( [ split[13], split[14], split[15] ] )
            X.append( split[16] )
        f.close()

        infile = files[i].replace("output", "input")

        X_xyz, X = np.asarray(X_xyz, dtype=np.float64), np.asarray(X, dtype=np.float64)
        #[ print(x) for x in X_xyz ]

        T = T*2.50

        inv_X = [ (1.0/x) for x in X ]

        avg_invX = np.average(inv_X)
        avg_X = np.average(X)
        Xm = getMolX( infile, avg_X )

        print("%6.4f   %10.8f   %10.8f" %(T, avg_X, avg_invX) )

        fit_x[i] = T
        fit_y[i] = avg_invX
        #fit_y[i] = avg_X

    A = np.vstack( [fit_x, np.ones(len(fit_x))] ).T
    coeffs = np.linalg.lstsq(A, fit_y, rcond=None)
    CW_pred = np.roots(coeffs[0])

    print("\n------------------------------------------------------------------------------\n")
    print("Fitted CW temperature = %10.8f (%4.3f)" %(abs(CW_pred), CW_pred) )
    print(coeffs)
    print("\n------------------------------------------------------------------------------\n")

    if plot:
       range_plot = [ -400.0, 400.0 ]
       plot_x = np.linspace( -400.0, 400, num=1000 )

       for i in range(len(plot_x)):
           print("%10.8f    %10.8f" %(plot_x[i], (coeffs[0]*plot_x[i])+coeffs[1] ) )

if __name__=="__main__":
   main()
