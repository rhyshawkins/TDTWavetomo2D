
import sys
import os
import argparse

import glob

import numpy
import matplotlib.pyplot as P

def load_residuals(filename):

    #
    # Residuals file has columns:
    #
    # 0 Station A
    # 1 Station B
    # 2 Distance km
    # 3 Raw Mean Residual
    # 4 Observered Velocity (to compare to raw)
    # 5 Normalised Mean Residual (should be N(0, 1) distributed)
    # 6 Mean Travel Time error seconds
    # 7 Observered Travel Time
    #
    

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    data = []
    
    for line in lines:
        t = line.split()
        
        A = t[0]
        B = t[1]
        nr = float(t[5])

        data.append((A, B, nr))

    return data

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file(s)')

    args = parser.parse_args()

    files = []
    if os.path.exists(args.input):
        files = [args.input]
    else:
        files = glob.glob(args.input + '-???')

    if len(files) == 0:
        print 'No files'
        sys.exit(-1)

    A, B, nr = zip(*load_residuals(files[0]))
    nr = numpy.array(nr)
    ndata = nr.size
    gnr = numpy.zeros((len(files), ndata))
    gnr[0, :] = nr
    
    for i, f in enumerate(files[1:]):

        _, _, nri = zip(*load_residuals(files[0]))

        gnr[i + 1, :] = numpy.array(nri)
        

    
    meannr = numpy.mean(gnr, axis = 0)

    fig, ax = P.subplots()

    ax.hist(meannr, bins = 50)


    indices = numpy.where(numpy.abs(meannr) > 3.0)[0]
    for i in indices:
        print A[i], B[i], meannr[i]

    print numpy.mean(meannr), numpy.std(meannr), ndata
    print '3.0: ', numpy.where(numpy.abs(meannr) > 3.0)[0].size, float(ndata) / 370.398
    print '3.3: ', numpy.where(numpy.abs(meannr) > 3.290527)[0].size, float(ndata) / 1000.0
    print '4.0: ', numpy.where(numpy.abs(meannr) > 4.0)[0].size, float(ndata) / 15787.0
    
    P.show()
        
