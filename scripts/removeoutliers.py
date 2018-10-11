
import sys
import os
import argparse

import glob

import numpy
import matplotlib.pyplot as P

def generate_paths(velocities):

    paths = {}

    for (sta, stb, _, _, _, _, _) in velocities:

        if sta not in paths:
            paths[sta] = []
        paths[sta].append(stb)

        if stb not in paths:
            paths[stb] = []
        paths[stb].append(sta)

    return paths

def eliminate_paths(paths):

    maxlen = -1
    maxkey = None
    for key, values in sorted(list(paths.items())):
        if len(values) > maxlen:
            maxlen = len(values)
            maxkey = key

    new_paths = {}
    for key, values in sorted(list(paths.items())):

        if key == maxkey:
            source = key
            destinations = values
        else:

            newvalues = []
            for v in values:
                if v != maxkey:
                    newvalues.append(v)

            if len(newvalues) > 0:
                new_paths[key] = newvalues

    return source, destinations, new_paths

def save_observations(filename, stations, freq, observations):

    paths = generate_paths(observations)
    traces = {}

    while len(paths) > 0:

        source, destinations, newpaths = eliminate_paths(paths)

        traces[source] = destinations
        
        paths = newpaths

    #
    # Output file
    #
    f = open(args.output, 'w')

    # Stations
    f.write('{:d}\n'.format(len(stations)))
    for key, (lon, lat) in sorted(list(stations.items())):
        f.write('{:s} {:15.9f} {:15.9f}\n'.format(key, lon, lat))

    # Frequencies
    f.write('1\n{:10.6f}\n'.format(freq))

    # Traces
    f.write('{:d}\n'.format(len(traces)))
    for key, pairs in sorted(list(traces.items())):
        f.write('{:s} {:d} {:s}\n'.format(key, len(pairs), ' '.join(pairs)))

    # Observations
    f.write('{:d}\n'.format(len(observations)))
    for (sta, stb, dist, vel, med, mod, err) in observations:
        f.write('{:s} {:s} {:15.9f}\n{:15.9f} {:15.9f} {:15.9f} {:15.9f}\n'.
          format(sta, stb, dist, vel, med, mod, err))

    f.close()
        

def load_observations(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    offset = 1
    nstations = int(lines[0])
    stations = {}
    for i in range(nstations):
        t = lines[offset].split()
        name = t[0]
        lon = float(t[1])
        lat = float(t[2])

        stations[name] = (lon, lat)
        offset = offset + 1

    nfreq = int(lines[offset])
    offset = offset + 1
    if nfreq != 1:
        raise Exception('No. frequencies must be 1 for this script.')

    freq = float(lines[offset])
    offset = offset + 1

    # Don't care about traces
    ntraces = int(lines[offset])
    offset = offset + 1 + ntraces

    nobs = int(lines[offset])
    offset = offset + 1

    data = []
    
    for i in range(nobs):

        t = lines[offset].split()
        A = t[0]
        B = t[1]
        distkm = float(t[2])
        
        mean, median, mode, err = map(float, lines[offset + 1].split())

        data.append((A, B, distkm, mean, median, mode, err))

        offset = offset + 2

    return freq, stations, data
        
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

def remove_observation(stations, observations, A, B):

    newobs = []
    nA = 0
    nB = 0
    for obs in observations:

        pA, pB, _, _, _, _, _ = obs

        if pA == A and pB == B:
            continue

        if pA == A or pB == A:
            nA = nA + 1

        if pB == B or pA == B:
            nB = nB + 1

        newobs.append(obs)

    newstations = {}
    for k, v in stations.items():
        if k == A and nA == 0:
            continue
        if k == B and nB == 0:
            continue

        newstations[k] = v
        
    return newstations, newobs

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file(s)')

    parser.add_argument('-d', '--data', type = str, required = True, help = 'Input data file')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output data file')

    parser.add_argument('-t', '--threshold', type = float, default = 3.0, help = 'Normed residual threshold')
    
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
    gnr = numpy.zeros((len(files), nr.size))
    gnr[0, :] = nr
    
    for i, f in enumerate(files[1:]):

        _, _, nri = zip(*load_residuals(files[0]))

        gnr[i + 1, :] = numpy.array(nri)
        

    
    meannr = numpy.mean(gnr, axis = 0)

    #
    # Load data
    #
    freq, stations, observations = load_observations(args.data)
    
    indices = numpy.where(numpy.abs(meannr) > args.threshold)[0]
    for i in indices:
        print 'Removing: ', A[i], B[i], meannr[i]

        stations, observations = remove_observation(stations, observations, A[i], B[i])

    print 'Removed : ', indices.size
    
    save_observations(args.output, stations, freq, observations)

