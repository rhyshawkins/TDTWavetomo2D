
import os
import glob
import argparse
import shutil

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file(s)')

    parser.add_argument('-p', '--path', type = str, required = True, help = 'Output path')
    parser.add_argument('-N', '--number', type = int, required = True, help = 'No. final models')

    args = parser.parse_args()

    files = glob.glob('%s-???' % args.input)

    if len(files) == 0:

        if os.path.exists(args.input):
            files = [args.input]

        else:

            print 'Error: no final model files found: %s' % args.input
            sys.exit(-1)

    M = len(files)
    for i in range(args.number):

        j = i % M

        src = files[j]
        dst = os.path.join(args.path, 'final_model.txt-%03d' % i)

        print dst
        shutil.copy(src, dst)

    
        

    
