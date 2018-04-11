//
//    TDTWavetomo2d : Software for the inversion of surface wave datasets using the
//    trans-dimensional tree approach using a wavelet parameterisation. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "volume.hpp"

static char short_options[] = "i:o:s:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"slice", required_argument, 0, 's'},
  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);
  
int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input;
  char *output;
  int slice;

  //
  // Defaults
  //
  input = nullptr;
  output = nullptr;
  slice = 0;

  //
  // Command line parameters
  //
  option_index = 0;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input = optarg;
      break;

    case 'o':
      output = optarg;
      break;

    case 's':
      slice = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (input == nullptr) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  double *volume;
  int width, height, depth;

  volume = volume_load(input, width, height, depth);
  if (volume == nullptr) {
    fprintf(stderr, "error: failed to read volume\n");
    return -1;
  }

  if (slice < 0 || slice >= depth) {
    fprintf(stderr, "error: slice %d out of range for depth %d\n", slice, depth);
    return -1;
  }
  
  if (volume_save_slice_as_text(output, volume, width, height, depth, slice) < 0) {
    fprintf(stderr, "error: failed to save slice\n");
    return -1;
  }

  delete [] volume;
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options are one or more of:\n"
	  "\n"
	  " -i|--input <filename>       Input volume filename\n"
	  " -o|--output <filename>      Output slice text filename\n"
	  "\n"
	  " -s|--slice <int>            Slice index\n"
	  "\n"
	  " -h|--help                   Show usage\n"
	  "\n",
	  pname);
}
