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

def main():

    plot = False

    try:
      #dirs = sorted([ ("%s%s"%(sys.argv[1], x)) for x in os.listdir( sys.argv[1] ) if ".py" not in x])
      dirs = sorted([ ("%s/%s"%(os.getcwd(), x)) for x in os.listdir( sys.argv[1] ) if ".py" not in x])
    except:
      print("Please enter directory!")
      sys.exit()

    T_scale = 0.0

    for i in range(len(dirs)):

        inner = sorted([ x for x in os.listdir( dirs[i] ) ])
        inner = [ ("%s/%s"%(dirs[i], x)) for x in inner ]

        files = [ ("%s/%s"%(x, "output")) for x in inner ]

        if "NiO" in dirs[i].split("/")[-1]:
            T_scale = 1.0
        else:
            T_scale = 2.5

        outfile_cur = ("%s.out"%((dirs[i].split("/")[-1])[3:]))
        write = []

        for j in range(len(files)):

            X_xyz, X = [], []
            T = np.float64(0.0)

            f=open(files[j], "r")
            [ next(f) for x in range(8) ]
            for k,line in enumerate(f):
                split = line.split()

                if k == 0:
                   T = np.float64(split[0] )

                X_xyz.append( [ split[13], split[14], split[15] ] )
                X.append( split[16] )
            f.close()

            X_xyz, X = np.asarray(X_xyz, dtype=np.float64), np.asarray(X, dtype=np.float64)
            T = T*2.50

            inv_X = [ (1.0/x) for x in X ]

            avg_invX = np.average(inv_X)
            avg_X = np.average(X)

            #print("%6.4f   %10.8f   %s" %(T, avg_X, "{:10.8E}".format(avg_invX)) )
            write.append( [ np.float64(T), "{:10.8E}".format(avg_X), "{:10.8E}".format(avg_invX) ] )

        f=open(outfile_cur, "w")
        f.write("T    X    X^-1\n")
        for l in range(len(write)):
            f.write("%4.3f\t %s\t %s\n" %(write[l][0], write[l][1], write[l][2])) 
        f.close()

    sys.exit()

    print("\n------------------------------------------------------------------------------\n")

if __name__=="__main__":
   main()
