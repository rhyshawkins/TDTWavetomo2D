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
#include <string.h>
#include <math.h>

#include <getopt.h>

extern "C" {
#include "chain_history.h"
  
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"

#include "slog.h"
  
};

#include "global.hpp"

struct user_data {
  int thincounter;
  int thin;
  int skip;
  
  int counter;
  
  int degree_max;
  int ncoeff;

  int **hist;
  int bins;
  double vmin;
  double vmax;

  wavetree3d_sub_t *wt;
  
};

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static int histogram_index(double v, double vmin, double vmax, int bins);

static char short_options[] = "d:l:i:o:t:s:b:z:Z:S:h";
static struct option long_options[] = {
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},

  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},

  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"bins", required_argument, 0, 'b'},
  {"vmin", required_argument, 0, 'z'},
  {"vmax", required_argument, 0, 'Z'},

  {"maxsteps", required_argument, 0, 'S'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  chain_history_t *ch;
  
  char *input_file;
  char *output_file;

  int degree_x;
  int degree_y;
  int degree_z;

  int thin;
  int skip;
  int maxsteps;

  int bins;
  double vmin;
  double vmax;

  FILE *fp_in;
  FILE *fp_out;

  struct user_data data;
  multiset_int_double_t *S_v;

  int i;
  int j;

  /*
   * Default values
   */
  fp_in = NULL;
  fp_out = NULL;
  degree_x = 4;
  degree_y = 4;
  degree_z = 4;
  
  input_file = NULL;
  output_file = NULL;
  
  bins = 1000;
  vmin = -3.0;
  vmax = 3.0;
  
  thin = 0;
  skip = 0;

  maxsteps = 1000000;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'x':
      degree_x = atoi(optarg);
      if (degree_x < 1) {
	fprintf(stderr, "error: invalid degree\n");
	return -1;
      }
      break;

    case 'y':
      degree_y = atoi(optarg);
      if (degree_y < 1) {
	fprintf(stderr, "error: invalid lateral degree\n");
	return -1;
      }
      break;

    case 'z':
      degree_z = atoi(optarg);
      if (degree_z < 1) {
	fprintf(stderr, "error: invalid lateral degree\n");
	return -1;
      }
      break;

    case 'i':
      input_file = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;
      
    case 'b':
      bins = atoi(optarg);
      if (bins < 1) {
	fprintf(stderr, "error: bins must be 1 or greater\n");
	return -1;
      }
      break;

    case 'S':
      maxsteps = atoi(optarg);
      if (maxsteps < 1000) {
	fprintf(stderr, "error: maxsteps should be 1000 or greater\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (input_file == NULL) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output_file == NULL) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  ch = chain_history_create(maxsteps);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }
  
  data.thincounter = 0;
  data.thin = thin;
  data.skip = skip;

  data.counter = 0;

  data.wt = wavetree3d_sub_create(degree_x, degree_y, degree_z, 0.0);
  if (data.wt == nullptr) {
    fprintf(stderr, "error: failed to create wavetree\n");
    return -1;
  }

  data.ncoeff = wavetree3d_sub_get_ncoeff(data.wt);
  data.degree_max = wavetree3d_sub_maxdepth(data.wt);

  data.bins = bins;
  data.vmin = vmin;
  data.vmax = vmax;
  data.hist =new int*[data.ncoeff];
  for (i = 0; i < data.ncoeff; i ++) {
    data.hist[i] = new int[data.bins];
    memset(data.hist[i], 0, sizeof(int) * data.bins);
  }
  
  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }
  
  fp_in = fopen(input_file, "r");
  if (fp_in == NULL) {
    fprintf(stderr, "error: failed to open input file\n");
    return -1;
  }

  /*
   * Process the chain history
   */
  while (!feof(fp_in)) {

    if (chain_history_read(ch,
			   (ch_read_t)fread,
			   fp_in) < 0) {
      if (feof(fp_in)) {
	break;
      }
      
      fprintf(stderr, "error: failed to read chain history\n");
      return -1;
    }

    if (chain_history_replay(ch,
			     S_v,
			     (chain_history_replay_function_t)process,
			     &data) < 0) {
      fprintf(stderr, "error: failed to replay\n");
      return -1;
    }
  }
  printf("%d records\n", data.counter);
  fclose(fp_in);

  fp_out = fopen(output_file, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "error: failed to open histogram file\n");
    return -1;
  }

  fprintf(fp_out, "%d %d\n", data.ncoeff, data.bins);
  fprintf(fp_out, "%.6f %.6f\n", data.vmin, data.vmax);
  
  for (j = 0; j < data.ncoeff; j ++) {
    for (i = 0; i < data.bins; i ++) {
      
      fprintf(fp_out, "%d ", data.hist[j][i]);
      
    }
    fprintf(fp_out, "\n");
  }
  
  fclose(fp_out);

  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  return 0;
}

static int process(int stepi,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;
  
  if ((d->thincounter >= d->skip) && (d->thin <= 1 || (d->thincounter % d->thin) == 0)) {

    for (int depth = 0; depth <= d->degree_max; depth ++) {

      int c = multiset_int_double_depth_count(S_v, depth);
      for (int i = 0; i < c; i ++) {
	int index;
	double value;
	
	if (multiset_int_double_nth_element(S_v, depth, i, &index, &value) < 0) {
	  ERROR("Failed to get nth element");
	  return -1;
	}

	int hi = histogram_index(value, d->vmin, d->vmax, d->bins);
	d->hist[index][hi] ++;
      }

    }
									     
  }
  d->thincounter ++;
  
  return 0;
}

static int histogram_index(double v, double vmin, double vmax, int bins)
{
  int i;
  
  i = (int)((double)bins * (v - vmin)/(vmax - vmin));

  if (i < 0) {
    return 0;
  }

  if (i > (bins - 1)) {
    return bins - 1;
  }

  return i;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -d|--degree-depth <int>    Number of layers as power of 2\n"
	  " -l|--degree-lateral <int>  Number of horizontal samples as power of 2\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  " -o|--output <file>               Output histograms for all coefficients\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -b|--bins <int>                  No. histogram bins\n"
	  " -z|--vmin <float>                Lower range for histogram\n"
	  " -Z|--vmax <float>                Upper range for histogram\n"
	  "\n"
	  " -S|--maxsteps <int>              Chain history max steps\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
