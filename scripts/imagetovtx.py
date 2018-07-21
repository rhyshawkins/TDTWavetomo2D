
import sys
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, default = None, help = 'Input image')
    parser.add_argument('-m', '--mean', type = float, default = None, help = 'Mean of image')

    parser.add_argument('-x', '--degreex', type = int, default = 5, help = 'x Degree of image')
    parser.add_argument('-y', '--degreey', type = int, default = 5, help = 'y Degree of image')

    parser.add_argument('-a', '--minlat', type = float, default = -1.0, help = 'Min lat')
    parser.add_argument('-A', '--maxlat', type = float, default = 1.0, help = 'Min lat')

    parser.add_argument('-n', '--minlon', type = float, default = -1.0, help = 'Min lon')
    parser.add_argument('-N', '--maxlon', type = float, default = 1.0, help = 'Min lon')

    parser.add_argument('-e', '--error', type = float, default = 0.1, help = 'Grid error')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')
    
    args = parser.parse_args()

    if args.mean is None and args.input is None:
        print 'You need to specify either a mean value or an image on the command line'
        sys.exit(-1)

    #
    # Construct image
    #
    width = 1 << args.degreex
    height = 1 << args.degreey
    
    if not args.input is None:
        image = numpy.loadtxt(args.mean)

        r, c = image.shape
        if r != height:
            height = r
            print 'Warning: height adjusted to %d' % height
        if c != width:
            width = c
            print 'Warning: width adjusted to %d' % width
    else:
        image = args.mean * numpy.ones((height, width))

    #
    # Write image in grid format
    #
    dlon = (args.maxlon - args.minlon)/float(width)
    gminlon = args.minlon - dlon/2.0

    dlat = (args.maxlat - args.minlat)/float(height)
    gmaxlat = args.maxlat + dlat/2.0
    
    f = open(args.output, 'w')
    f.write('%d %d\n' % (height, width))
    f.write('%10.6f %10.6f\n' % (gmaxlat, gminlon))
    f.write('%15.9f %16.9f\n' % (dlat, dlon))

    for i in range(width + 2):

        ii = i - 1
        if (ii < 0):
            ii = 0
        elif (ii == width):
            ii = width - 1
            
        for j in range(height + 2):

            jj = j - 1
            if (jj < 0):
                jj = 0
            elif (jj == height):
                jj = height - 1
                
            
            f.write('%10.6f %10.6f\n' % (image[height - jj - 1, ii], args.error))

    f.close()
    
    
    
    
