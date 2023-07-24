import os 
import sys
import numpy as np
import subprocess
import glob

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

def getdata(in_file, size):
    T = np.zeros( size, dtype=np.float64 )
    M_a, M_b = np.zeros( size, dtype=np.float64 ), np.zeros( size, dtype=np.float64 )
    X = np.zeros( size, dtype=np.float64 )
    f=open(in_file, "r")
    [ next(f) for x in range(8) ]
    for i,line in enumerate(f):
        split = line.split()
        T[i] = np.float64(split[0])
        M_a[i] = np.float64(split[4])
        M_b[i] = np.float64(split[8])
        X[i] = np.float64(split[16])
    f.close()
    return T, M_a, M_b, X

def arrsize(in_file):
    size = 0
    f=open(in_file, "r")
    [ next(f) for x in range(8) ]
    for line in f:
        size +=1
    f.close()
    return size

def GausSmooth( T, X, sigma, Tst, numst ):
    Tst, Tend = np.float64(Tst), np.float64(0.0)
    dT = np.float64( (Tend-Tst)/numst )
    sigma = np.float64(sigma)
    den, top = np.float64(0.0), np.float64(0.0)
    X_filter = np.zeros( numst+1, dtype=np.float64 )
    T_new = np.zeros( numst+1, dtype=np.float64 )
    for i in range( numst +1):
        if i:
          T_cur = T_cur - dT
        else:
          T_cur = Tend
        den = np.float64(0.0)
        top = np.float64(0.0)
        for j in range(len(T)):
            den += np.exp( - np.power( (T_cur - T[j]), 2)/np.power(sigma,2) )
            top += X[j]*np.exp( - np.power( (T_cur - T[j]), 2)/np.power(sigma,2) )
        try:           
           X_filter[i] = top/den
        except:
           X_filter[i] = 0.0
        T_new[i] = T_cur
    return T_new, X_filter

def main():

    smooth = False
    sig = 2.5
    fmt = "{:8.6E}"

    try: 
       in_dir = sys.argv[1]
    except:
       print("\n Enter a directory! \n")
       sys.exit()

    dirs = sorted([ ("%s/%s/"%(os.getcwd(), x)) for x in os.listdir( in_dir ) if ".py" not in x ])
    Tst, numst = 0.0, 0.0

    for i in range(len( dirs )):
        dump = sorted( [ ("%s%s" %(dirs[i], x) ) for x in os.listdir( dirs[i] ) ] )
        dummy = pyGrep( [ "atomic-spin-moment", glob.glob( ("%s/*.mat"%(dirs[i]) ) )[0] ] )
        T_scale = np.float64( dummy.split()[0].split("=")[-1] )
        out_cur = ("%s/output"%(dirs[i]))
        size = arrsize(out_cur)

        T, M_a, M_b, X = getdata(out_cur, size)

        if "NiO" in dirs[i].split("/")[-2]:
           f1_cur = ("%s_raw.out"%( (dirs[i].split("/")[-2])[3:] ))
           f2_cur = ("%s_smooth.out"%( (dirs[i].split("/")[-2])[3:] ))
           Tst, numst = 1000.0, 2000
            
        else:
           f1_cur = ("%s_%s_raw.out" %(dirs[i].split("_")[-2].split("/")[0], dirs[i].split("_")[-1])).replace("/","")
           f2_cur = ("%s_%s_smooth.out" %(dirs[i].split("_")[-2].split("/")[0], dirs[i].split("_")[-1])).replace("/","")
           Tst, numst = 200.0, 400

        T_new, X_new = GausSmooth( T, X, sig, Tst, numst )

        f1=open(f1_cur, "w")
        for k in range(size):
            M_avg = np.float64( ( M_a[k] + M_b[k] )/2.0 )
            f1.write("%4.3f\t%s\t%s\n" %(T[k]*T_scale, fmt.format(M_avg), "{:14.12E}".format(1.0/X[k])) )
        f1.close()

        f2=open(f2_cur, "w")
        for l in range(len( X_new )):
            f2.write("%4.3f\t%s\n"%(T_new[l]*T_scale, "{:14.12E}".format(1.0/X_new[l]) ))
        f2.close()
            
if __name__ == "__main__":
   main()
