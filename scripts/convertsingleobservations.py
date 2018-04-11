#!/usr/bin/env python3
#-*- coding:utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

import argparse

def load_stations(filename):

    f = open(filename, 'r')
    lines = f.readlines()

    stations = {}
    for line in lines:
        t = line.split()

        name = t[0]
        lon = float(t[1])
        lat = float(t[2])

        stations[name] = (lon, lat)

    return stations

def load_velocities(filename, stations):
    
    f = open(filename, 'r')
    lines = f.readlines()

    velocities = []

    for line in lines:
        t = line.split()

        sta = t[0]
        stb = t[1]
        dist = float(t[2])
        vel = float(t[3])
        err = float(t[4])

        if sta in stations:
            if stb in stations:

                velocities.append((sta, stb, dist, vel, err))
            else:
                print('warning: unknown station {:s} in velocities'.format(stb))
                
        else:
            print('warning: unknown station {:s} in velocities'.format(sta))

    return velocities

def generate_paths(velocities):

    paths = {}

    for (sta, stb, _, _, _) in velocities:

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
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--stations', type = str, required = True, help = 'Station list')
    parser.add_argument('-v', '--velocities', type = str, required = True, help = 'Interstation velocities')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('-f', '--frequency', type = float, default = 1.0, help = 'Frequency (Hz)')
    
    args = parser.parse_args()

    #
    # Load stations
    #
    stations = load_stations(args.stations)

    velocities = load_velocities(args.velocities, stations)

    #
    # Generate traces
    #
    paths = generate_paths(velocities)
    
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
    f.write('1\n{:10.6f}\n'.format(args.frequency))

    # Traces
    f.write('{:d}\n'.format(len(traces)))
    for key, pairs in sorted(list(traces.items())):
        f.write('{:s} {:d} {:s}\n'.format(key, len(pairs), ' '.join(pairs)))

    # Observations
    f.write('{:d}\n'.format(len(velocities)))
    for (sta, stb, dist, vel, err) in velocities:
        f.write('{:s} {:s} {:15.9f}\n{:15.9f} {:15.9f} {:15.9f} {:15.9f}\n'.
          format(sta, stb, dist, vel, vel, vel, err))

    f.close()
