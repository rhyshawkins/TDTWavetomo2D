

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    args = parser.parse_args()

    f = open(args.input, 'r')
    lines = f.readlines()
    f.close()
    
    o = 0

    nstations = int(lines[o])
    o = o + nstations + 1
    
    f = open(args.output, 'w')
    for i in range(nstations):

        f.write(lines[i + 1])

    f.close()
    

    
