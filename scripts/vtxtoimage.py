
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input vtx file')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output image file')

    args = parser.parse_args()

    f = open(args.input, 'r')
    lines = f.readlines()
    f.close()
    
    height, width = map(int, lines[0].split())
    latmin, lonmin = map(float, lines[1].split())
    dlat, dlon = map(float, lines[2].split())

    image = numpy.zeros((height - 2, width - 2))
    error = numpy.zeros((height - 2, width - 2))
    
    ii = 0
    jj = 0
    
    for line in lines[3:]:

        if line.strip() == '':
            continue

        v, err = map(float, line.split())

        if (jj > 0 and ii > 0 and jj < (height - 1) and ii < (width - 1)):
            image[jj - 1, ii - 1] = v
            error[jj - 1, ii - 1] = err

        jj = jj + 1
        if jj == height:
            jj = 0
            ii = ii + 1

    numpy.savetxt(args.output, image)
    numpy.savetxt(args.output + '.error', error)
    

