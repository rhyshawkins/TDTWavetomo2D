
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('-s', '--slice', type = int, default = 0, help = 'Frequency slice')
    
    args = parser.parse_args()

    f = open(args.input, 'r')
    lines = f.readlines()
    o = 0

    nstations = int(lines[o])
    o = o + nstations + 1
    print nstations
    
    nfreq = int(lines[o])
    o = o + nfreq + 1
    print nfreq

    ntraces = int(lines[o])
    o = o + ntraces + 1
    print ntraces

    nobs = int(lines[o])
    o = o + 1
    print nobs

    velocities = []
    for i in range(nobs):
        t = lines[o].split()

        sta = t[0]
        stb = t[1]
        dist = float(t[2])
        o = o + 1

        di = o + args.slice
        v, _, _, err = map(float, lines[di].split())
        
        o = o + nfreq

        velocities.append((sta, stb, dist, v, err))

    f.close()

    f = open(args.output, 'w')
    for (sta, stb, dist, v, err) in velocities:

        f.write('%s %s %15.9f %15.9f %15.9f\n' % (sta, stb, dist, v, err))

    f.close()
    

    
