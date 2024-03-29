
import sys
import os
import argparse

import glob

import numpy
import matplotlib.pyplot as P

from mpl_toolkits.basemap import Basemap

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

def filter_files(files, chains):
    if chains <= 0:
        return files

    else:
        return filter(lambda x: int(x[-3:]) < chains, files)

def load_stations(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    N = int(lines[0])

    stations = {}
    lons = []
    lats = []
    for i in range(N):
        t = lines[i + 1].split()

        name = t[0]
        lon = float(t[1])
        lat = float(t[2])

        stations[name] = (lon, lat)

        lons.append(lon)
        lats.append(lat)

    

    return stations, min(lons), max(lons), min(lats), max(lats)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file(s)')

    parser.add_argument('-c', '--chains', type = int, default = -1, help = 'No. chains (at T = 1)')

    parser.add_argument('-t', '--threshold', type = float, default = 3.0, help = 'Outlier threshold')

    parser.add_argument('-d', '--data', type = str, default = None, help = 'Plot stations')

    args = parser.parse_args()

    files = []
    if os.path.exists(args.input):
        files = [args.input]
    else:
        files = filter_files(glob.glob(args.input + '-???'), args.chains)
        
        

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


    indices = numpy.where(numpy.abs(meannr) > args.threshold)[0]
    if len(indices) > 0:
        t = zip(numpy.abs(meannr[indices]), indices)
        t.sort()
        _, sindices = zip(*t)

        print 'Outliers'
        print '-' * 60
        for i in sindices:
            print '%20s %20s %10.6f' % (A[i], B[i], meannr[i])
        print '-' * 60

    print 'Mean Normed Residual: %16.9e' % numpy.mean(meannr)
    print 'Std. Dev            : %16.9e (%d)' % (numpy.std(meannr), ndata)

    print '%19s : %10s : %10s' % ('Threshold', 'Count', 'Expected')
    print '%19.3f : %10d : %10.6f' % (3.0,
                                      numpy.where(numpy.abs(meannr) > 3.0)[0].size,
                                      float(ndata) / 370.398)
    
    print '%19.3f : %10d : %10.6f' % (3.291,
                                      numpy.where(numpy.abs(meannr) > 3.290527)[0].size,
                                      float(ndata) / 1000.0)

    print '%19.3f : %10d : %10.6f' % (4.0,
                                      numpy.where(numpy.abs(meannr) > 4.0)[0].size,
                                      float(ndata) / 15787.0)

    if len(indices) > 0 and (not (args.data is None)):

        stations, minlon, maxlon, minlat, maxlat = load_stations(args.data)

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

            if meannr[i] < 0.0:
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
        
