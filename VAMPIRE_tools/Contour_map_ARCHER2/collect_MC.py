import os
import sys
import numpy as np

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

def main():

    in_file = sys.argv[1]
    size = arrsize(in_file)
    T, M_a, M_b, X = getdata(in_file, size)
    T_scale = 2.5

    for i in range(size):
        #M_avg = np.float64( ( M_a[i] + M_b[i] )/2.0 )
        print("%4.3f    %10.8f    %10.8f" %(T[i]*2.5, M_avg, X[i]) )

if __name__ == "__main__":
   main()
