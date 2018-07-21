
import argparse

def load_data(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    nstations = int(lines[0])
    o = 1
    stations = {}
    for i in range(nstations):
        t = lines[o].split()

        lon = float(t[1])
        lat = float(t[2])

        stations[t[0]] = (lon, lat)

        o = o + 1

    nfrequencies = int(lines[o])
    o = o + 1
    frequencies = []
    for i in range(nfrequencies):
        frequencies.append(float(lines[o]))
        o = o + 1

    ntraces = int(lines[o])
    o = o + 1
    traces = {}
    for i in range(ntraces):
        t = lines[o].split()

        traces[t[0]] = t[2:]
        o = o + 1

    nobservations = int(lines[o])
    o = o + 1
    observations = {}
    for i in range(nobservations):
        t = lines[o].split()
        A = t[0]
        B = t[1]
        distkm = float(t[2])
        o = o + 1

        data = []
        for j in range(nfrequencies):
            mean, med, mode, err = map(float, lines[o].split())
            o = o + 1
            data.append((mean, med, mode, err))

        if (traces.has_key(A) and B in traces[A]):
            if not observations.has_key(A):
                observations[A] = {}

            observations[A][B] = (distkm, data)
        else:
            if not observations.has_key(B):
                observations[B] = {}

            observations[B][A] = (distkm, data)

    return stations, frequencies, traces, observations

def write_stations(stations, order, filename):

    f = open(filename, 'w')

    f.write('%d\n' % len(order))
    for nm in order:

        lon, lat = stations[nm]
        f.write('%10.6f %10.6f\n' % (lat, lon))

    f.close()

def write_station_order(stations, order, filename):
    
    f = open(filename, 'w')

    f.write('%d\n' % len(order))
    for nm in order:

        lon, lat = stations[nm]
        f.write('%s %10.6f %10.6f\n' % (nm, lat, lon))

    f.close()


def write_observations(observations, order, slicei, filename):

    f = open(filename, 'w')
    
    for sname in order:

        haspaths = observations.has_key(sname)
        
        for rname in order:

            if haspaths and observations[sname].has_key(rname):

                distkm, data = observations[sname][rname]
                mean, _, _, err = data[slicei]

                ttime = distkm/mean

                #
                # Approximate conversion of error in velocity to travel time
                #
                dt1 = distkm/(mean - err) - ttime
                dt2 = ttime - distkm/(mean + err)
                terr = (dt1 + dt2)/2.0

                f.write('1 %10.6f %10.6f\n' % (ttime, terr))

            else:
                f.write('0 0.0 0.0\n')

    f.close()

                        
                

                
            
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input observations')
    parser.add_argument('-o', '--output-prefix', type = str, required = True, help = 'Base name for output files')

    parser.add_argument('-s', '--slice', type = int, default = 0, help = 'Frequency slice')

    args = parser.parse_args()

    stations, frequencies, traces, observations = load_data(args.input)

    #
    # First write the sources and receivers (same file), we generate a consistent
    # station order here
    #
    order = stations.keys()
    order.sort()
    
    write_stations(stations, order, '%s-sources.dat' % args.output_prefix)
    write_stations(stations, order, '%s-receivers.dat' % args.output_prefix)
    write_station_order(stations, order, '%s-order.dat' % args.output_prefix)

    #
    # Write the observations
    #
    write_observations(observations, order, args.slice, '%s-otimes.dat' % args.output_prefix)
    
