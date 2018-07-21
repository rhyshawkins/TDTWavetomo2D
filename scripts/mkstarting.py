
import os
import shutil
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output directory')

    parser.add_argument('-n', '--np', type = int, default = 8, help = 'No. input files')

    args = parser.parse_args()

    try:
        os.stat(args.output)

    except OSError, e:
        os.mkdir(args.output)

    for i in range(args.np):

        filename = os.path.join(args.output, 'final-model.txt-%03d' % i)

        shutil.copy(args.input, filename)

        
    
    
