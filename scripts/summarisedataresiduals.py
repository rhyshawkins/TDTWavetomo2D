
import sys
import os
import argparse

import glob

import numpy
import matplotlib.pyplot as P

from mpl_toolkits.basemap import Basemap

def filter_files(files, chains):
    if chains <= 0:
        return files

    else:
        return filter(lambda x: int(x[-3:]) < chains, files)

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

def stationrange(stations):

    minlon = 180.0
    maxlon = -180.0
    minlat = 90.0
    maxlat = -90.0

    for _, (lon, lat) in stations.items():

        minlon = min([minlon, lon])
        maxlon = max([maxlon, lon])
        minlat = min([minlat, lat])
        maxlat = max([maxlat, lat])

    return minlon, maxlon, minlat, maxlat

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file(s)')

    parser.add_argument('-t', '--threshold', type = float, default = 3.0, help = 'Outlier threshold')

    args = parser.parse_args()

    _, stations, data = load_observations(args.input)

    A, B, distkm, mean, _, _, err = zip(*data)

    vmean = numpy.mean(mean)

    res = (numpy.array(mean) - vmean)/numpy.array(err)
    
    fig, ax = P.subplots()

    ax.hist(res, bins = 50)


    indices = numpy.where(numpy.abs(res) > args.threshold)[0]
    if len(indices) > 0:
        t = zip(numpy.abs(res[indices]), indices)
        t.sort()
        _, sindices = zip(*t)

        print 'Outliers'
        print '-' * 60
        for i in sindices:
            print '%20s %20s %10.6f' % (A[i], B[i], res[i])
        print '-' * 60

    ndata = len(A)
    print 'Mean Normed Residual: %16.9e' % numpy.mean(res)
    print 'Std. Dev            : %16.9e (%d)' % (numpy.std(res), ndata)

    print '%19s : %10s : %10s' % ('Threshold', 'Count', 'Expected')
    print '%19.3f : %10d : %10.6f' % (3.0,
                                      numpy.where(numpy.abs(res) > 3.0)[0].size,
                                      float(ndata) / 370.398)
    
    print '%19.3f : %10d : %10.6f' % (3.291,
                                      numpy.where(numpy.abs(res) > 3.290527)[0].size,
                                      float(ndata) / 1000.0)

    print '%19.3f : %10d : %10.6f' % (4.0,
                                      numpy.where(numpy.abs(res) > 4.0)[0].size,
                                      float(ndata) / 15787.0)

    if len(indices) > 0:

        minlon, maxlon, minlat, maxlat = stationrange(stations)
        
        dlon = 0.10 * (maxlon - minlon)
        minlon = minlon - dlon
        maxlon = maxlon + dlon
        
        dlat = 0.10 * (maxlat - minlat)
        minlat = minlat - dlat
        maxlat = maxlat + dlat

        fig, ax = P.subplots()

        m = Basemap(llcrnrlon = minlon,
                    llcrnrlat = minlat,
                    urcrnrlon = maxlon,
                    urcrnrlat = maxlat,
                    rsphere = (6378137.0, 6356752.3142),
                    resolution = 'i',
                    projection = 'merc',
                    ax = ax)

        m.fillcontinents(color = 'coral', lake_color = 'aqua')
        m.drawcoastlines(linewidth = 0.25)
        m.drawmapboundary(fill_color = 'aqua')
        
        unique_stations = {}
        
        for i in indices:

            lon0, lat0 = stations[A[i]]
            lon1, lat1 = stations[B[i]]

            unique_stations[A[i]] = (lon0, lat0)
            unique_stations[B[i]] = (lon1, lat1)

            if res[i] < 0.0:
                m.drawgreatcircle(lon0, lat0, lon1, lat1, color = 'red', zorder = 50)
            else:
                m.drawgreatcircle(lon0, lat0, lon1, lat1, color = 'blue', zorder = 50)

        for name, (lon, lat) in unique_stations.items():

            x, y = m(lon, lat)
            ax.scatter([x], [y], marker = '*', color = (0, 0, 0, 1), zorder = 100)
            ax.text(x, y - 1.5e4, name,
                    verticalalignment = 'top',
                    horizontalalignment = 'center')
    
    P.show()
        
