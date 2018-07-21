
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

    image = numpy.zeros((height, width))
    #error = numpy.zeros((height, width))
    
    i = 0
    
    for line in lines[3:]:

        if line.strip() == '':
            continue

        v = float(line)

        jj = i % (height + 2)
        ii = i / (height + 2)
        i = i + 1
        
        if (jj > 0 and ii > 0 and jj < (height + 1) and ii < (width + 1)):

            j = jj - 1
            image[height - j - 1, ii - 1] = v
            
            #error[jj - 1, ii - 1] = err

    numpy.savetxt(args.output, image)
    #numpy.savetxt(args.output + '.error', error)
    

