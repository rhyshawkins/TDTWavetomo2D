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
  
};

#include "global.hpp"
#include "wavetomo2dexception.hpp"
#include "wavetomo2dutil.hpp"
#include "volume.hpp"

struct user_data {
  int thincounter;
  int thin;
  int skip;
  
  int counter;

  int depthmax;
  int width;
  int height;
  int depth;
  int slicestride;
  int size;

  double *zoffset;
  int zoffset_n;

  double *mean;
  double *variance;

  int **hist;
  int bins;
  double vmin;
  double vmax;

  double *model;
  double *workspace;

  generic_lift_inverse1d_step_t hwaveletf;
  generic_lift_inverse1d_step_t vwaveletf;

  wavetree3d_sub_t *wt;

};

static const double CREDIBLE_INTERVAL = 0.95;

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static int histogram_index(double v, double vmin, double vmax, int bins);

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins);
static double median_from_histogram(int *hist, double vmin, double vmax, int bins);
static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double hpd_from_histogram(int *hist, double vmin, double vmax, int bins, double hpd_interval, double &hpd_min, double &hpd_max);

static char short_options[] = "x:y:z:i:Z:o:D:t:s:m:M:c:C:g:p:P:Q:b:v:V:S:w:W:F:h";
static struct option long_options[] = {
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},

  {"input", required_argument, 0, 'i'},
  {"zoffset", required_argument, 0, 'Z'},
  {"output", required_argument, 0, 'o'},
  
  {"stddev", required_argument, 0, 'D'},
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"mode", required_argument, 0, 'm'},
  {"median", required_argument, 0, 'M'},
  {"credible-min", required_argument, 0, 'c'},
  {"credible-max", required_argument, 0, 'C'},
  {"histogram", required_argument, 0, 'g'},

  {"hpd-min", required_argument, 0, 'p'},
  {"hpd-max", required_argument, 0, 'P'},
  {"hpd-range", required_argument, 0, 'Q'},

  {"bins", required_argument, 0, 'b'},
  {"vmin", required_argument, 0, 'v'},
  {"vmax", required_argument, 0, 'V'},

  {"maxsteps", required_argument, 0, 'S'},

  {"wavelet-lateral", required_argument, 0, 'w'},
  {"wavelet-depth", required_argument, 0, 'W'},

  {"filter-depth", required_argument, 0, 'F'},

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
  char *zoffset_file;
  char *output_file;
  char *stddev_file;

  int degree_x;
  int degree_y;
  int degree_z;

  int thin;
  int skip;
  int maxsteps;

  char *mode_file;
  char *median_file;
  char *credible_min;
  char *credible_max;
  char *histogram;
  char *hpd_min;
  char *hpd_max;
  char *hpd_range;

  int bins;
  double vmin;
  double vmax;

  FILE *fp_in;
  FILE *fp_out;

  struct user_data data;
  multiset_int_double_t *S_v;

  int i;
  int j;

  int credible_drop;

  int waveleth;
  int waveletv;
  int depthmax;

  /*
   * Default values
   */
  fp_in = NULL;
  fp_out = NULL;
  
  degree_x = 4;
  degree_y = 4;
  degree_z = 4;
  
  input_file = NULL;
  zoffset_file = NULL;
  output_file = NULL;
  stddev_file = NULL;
  
  mode_file = NULL;
  median_file = NULL;
  credible_min = NULL;
  credible_max = NULL;
  histogram = NULL;

  hpd_min = NULL;
  hpd_max = NULL;
  hpd_range = NULL;
  
  bins = 1000;
  vmin = 2.0;
  vmax = 4.0;
  
  thin = 0;
  skip = 0;

  maxsteps = 1000000;

  waveleth = 0;
  waveletv = 0;

  depthmax = -1;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'x':
      degree_x = atoi(optarg);
      if (degree_x < 1) {
	fprintf(stderr, "error: invalid x degree\n");
	return -1;
      }
      break;

    case 'y':
      degree_y = atoi(optarg);
      if (degree_y < 1) {
	fprintf(stderr, "error: invalid y degree\n");
	return -1;
      }
      break;

    case 'z':
      degree_z = atoi(optarg);
      if (degree_z < 1) {
	fprintf(stderr, "error: invalid z degree\n");
	return -1;
      }
      break;

    case 'i':
      input_file = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'D':
      stddev_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;
      
    case 'm':
      mode_file = optarg;
      break;
      
    case 'M':
      median_file = optarg;
      break;

    case 'c':
      credible_min = optarg;
      break;

    case 'C':
      credible_max = optarg;
      break;

    case 'g':
      histogram = optarg;
      break;

    case 'p':
      hpd_min = optarg;
      break;

    case 'P':
      hpd_max = optarg;
      break;

    case 'Q':
      hpd_range = optarg;
      break;

    case 'b':
      bins = atoi(optarg);
      if (bins < 1) {
	fprintf(stderr, "error: bins must be 1 or greater\n");
	return -1;
      }
      break;

    case 'v':
      vmin = atof(optarg);
      break;

    case 'V':
      vmax = atof(optarg);
      break;

    case 'S':
      maxsteps = atoi(optarg);
      if (maxsteps < 1000) {
	fprintf(stderr, "error: maxsteps should be 1000 or greater\n");
	return -1;
      }
      break;

    case 'w':
      waveletv = atoi(optarg);
      if (waveletv < 0 || waveletv > Global::WAVELET_MAX) {
	fprintf(stderr, "error: vertical wavelet must be between 0 and %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'W':
      waveleth = atoi(optarg);
      if (waveleth < 0 || waveleth > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be between 0 and %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'Z':
      zoffset_file = optarg;
      break;

    case 'F':
      depthmax = atoi(optarg);
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

  data.depthmax = depthmax;
  data.width = wavetree3d_sub_get_width(data.wt);
  data.height = wavetree3d_sub_get_height(data.wt);
  data.depth = wavetree3d_sub_get_depth(data.wt);
  data.slicestride = data.width * data.height;
  data.size = wavetree3d_sub_get_size(data.wt);

  printf("Image: %d x %d x %d (%d)\n", data.width, data.height, data.depth, data.size);

  if (zoffset_file != NULL) {
    data.zoffset = load_zoffset(zoffset_file, data.zoffset_n);
    if (data.zoffset == NULL) {
      fprintf(stderr, "error: failed to load zoffset file\n");
      return -1;
    }

    if (data.zoffset_n != data.depth) {
      fprintf(stderr, "error: invalid zoffset file (%d != %d)\n", data.zoffset_n, data.depth);
      return -1;
    }
    
  } else {
    data.zoffset = NULL;
    data.zoffset_n = 0;
  }
  
  data.mean = new double[data.size];
  memset(data.mean, 0, sizeof(double) * data.size);

  data.bins = bins;
  data.vmin = vmin;
  data.vmax = vmax;
  data.hist =new int*[data.size];
  for (i = 0; i < data.size; i ++) {
    data.hist[i] = new int[data.bins];
    memset(data.hist[i], 0, sizeof(int) * data.bins);
  }
  
  data.variance = new double[data.size];
  memset(data.variance, 0, sizeof(double) * data.size);
  
  data.model = new double[data.size];

  int workspacesize = data.width;
  if (data.height > workspacesize) {
    workspacesize = data.height;
  }
  if (data.depth > workspacesize) {
    workspacesize = data.depth;
  }
  data.workspace = new double[workspacesize];

  data.vwaveletf = Global::wavelet_inverse_function_from_id(waveletv);
  data.hwaveletf = Global::wavelet_inverse_function_from_id(waveleth);

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

  /*
   * Mean output
   */
  if (volume_save(output_file, data.mean, data.width, data.height, data.depth) < 0) {
    fprintf(stderr, "error: failed to save mean\n");
    return -1;
  }

  for (i = 0; i < data.size; i ++) {
    data.variance[i] = sqrt(data.variance[i]/(double)(data.counter - 1));
  }

  /*
   * Std. Deviation output
   */
  if (stddev_file != NULL) {
    if (volume_save(output_file, data.variance, data.width, data.height, data.depth) < 0) {
      fprintf(stderr, "error: failed to save std dev\n");
      return -1;
    }
  }

  /*
   * Mode output
   */
  if (mode_file != NULL) {
    //
    // Construct mode from histogram and store in mean
    //
    for (i = 0; i < data.size; i ++) {
      data.mean[i] = mode_from_histogram(data.hist[i],
					 data.vmin, data.vmax, data.bins);
    }

    if (volume_save(mode_file, data.mean, data.width, data.height, data.depth) < 0) {
      fprintf(stderr, "error: failed to save mode\n");
      return -1;
    }
  }
  
  /*
   * Median output
   */
  if (median_file != NULL) {
    fp_out = fopen(median_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open median file\n");
      return -1;
    }

    for (j = 0; j < data.height; j ++) {
      for (i = 0; i < data.width; i ++) {
	fprintf(fp_out, "%10.6f ", median_from_histogram(data.hist[j*data.width + i],
							 data.vmin, data.vmax, data.bins));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  /*
   * Credible Min
   */
  credible_drop = (int)(((double)data.counter * (1.0 - CREDIBLE_INTERVAL))/2.0);
  
  if (credible_min != NULL) {
    fp_out = fopen(credible_min, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open credible min file\n");
      return -1;
    }

    for (j = 0; j < data.height; j ++) {
      for (i = 0; i < data.width; i ++) {
	fprintf(fp_out, "%10.6f ", head_from_histogram(data.hist[j*data.width + i],
						       data.vmin, data.vmax, data.bins,
						       credible_drop));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }
  
  /*
   * Credible Max
   */
  if (credible_max != NULL) {
    fp_out = fopen(credible_max, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open credible max file\n");
      return -1;
    }

    for (j = 0; j < data.height; j ++) {
      for (i = 0; i < data.width; i ++) {
	fprintf(fp_out, "%10.6f ", tail_from_histogram(data.hist[j*data.width + i],
						       data.vmin, data.vmax, data.bins,
						       credible_drop));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  if (histogram != NULL) {
    fp_out = fopen(histogram, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open histogram file\n");
      return -1;
    }

    fprintf(fp_out, "%d %d\n", data.size, data.bins);
    fprintf(fp_out, "%.6f %.6f\n", data.vmin, data.vmax);

    for (j = 0; j < data.size; j ++) {
      for (i = 0; i < data.bins; i ++) {

	fprintf(fp_out, "%d ", data.hist[j][i]);

      }
      fprintf(fp_out, "\n");
    }
    
    fclose(fp_out);
  }

  if (hpd_range != NULL || hpd_min != NULL || hpd_max != NULL) {

    FILE *fpr = NULL;
    FILE *fpmin = NULL;
    FILE *fpmax = NULL;

    if (hpd_range != NULL) {
      fpr = fopen(hpd_range, "w");
      if (fpr == NULL) {
	fprintf(stderr, "error: failed to open hpd range file\n");
	return -1;
      }
    }
    
    if (hpd_min != NULL) {
      fpmin = fopen(hpd_min, "w");
      if (fpmin == NULL) {
	fprintf(stderr, "error: failed to open hpd min file\n");
	return -1;
      }
    }

    if (hpd_max != NULL) {
      fpmax = fopen(hpd_max, "w");
      if (fpmax == NULL) {
	fprintf(stderr, "error: failed to open hpd max file\n");
	return -1;
      }
    }

    for (j = 0; j < data.height; j ++) {
      for (i = 0; i < data.width; i ++) {
	double hmin, hmax;

	double hrange = hpd_from_histogram(data.hist[j*data.width + i],
					   data.vmin, data.vmax, data.bins,
					   CREDIBLE_INTERVAL,
					   hmin, hmax);

	if (fpr != NULL) {
	  fprintf(fpr, "%10.6f ", hrange);
	}

	if (fpmin != NULL) {
	  fprintf(fpmin, "%10.6f ", hmin);
	}
	  
	if (fpmax != NULL) {
	  fprintf(fpmax, "%10.6f ", hmax);
	}
      }

      if (fpr != NULL) {
	fprintf(fpr, "\n");
      }

      if (fpmin != NULL) {
	fprintf(fpmin, "\n");
      }

      if (fpmax != NULL) {
	fprintf(fpmax, "\n");
      }
    }

    if (fpr != NULL) {
      fclose(fpr);
    }
    
    if (fpmin != NULL) {
      fclose(fpmin);
    }

    if (fpmax != NULL) {
      fclose(fpmax);
    }
  }    

  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  delete [] data.mean;
  delete [] data.variance;
  delete [] data.model;
  delete [] data.workspace;
  
  return 0;
}

static int process(int stepi,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;
  double delta;
  int i;
  int j;
  int hi;
  
  if ((d->thincounter >= d->skip) && (d->thin <= 1 || (d->thincounter % d->thin) == 0)) {
    
    memset(d->model, 0, sizeof(double) * d->size);

    if (d->depthmax < 0) {
      if (wavetree3d_sub_set_from_S_v(d->wt, S_v) < 0) {
	fprintf(stderr, "process: failed to set wavetree (sub)\n");
	return -1;
      }
    } else {
      if (wavetree3d_sub_set_from_S_v_filtered(d->wt, S_v, d->depthmax) < 0) {
	fprintf(stderr, "process: failed to set wavetree filtered (sub)\n");
	return -1;
      }
    }
    
    if (wavetree3d_sub_map_to_array(d->wt, d->model, d->size) < 0) {
      fprintf(stderr, "process: failed to map to array\n");
      return -1;
    }

    
    if (generic_lift_inverse3d(d->model,
			       d->width,
			       d->height,
			       d->depth,
			       d->width,
			       d->slicestride,
			       d->workspace,
			       d->hwaveletf,
			       d->hwaveletf,
			       d->vwaveletf,
			       1) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    //
    // Added zoffset if available
    //
    if (d->zoffset != NULL) {
      for (j = 0; j < d->depth; j ++) {

	for (i = 0; i < d->slicestride; i ++) {

	  d->model[j*d->slicestride + i] += d->zoffset[j];

	}
      }
    }
    
    /*
     * Update mean/variance calculation
     */
    d->counter ++;
    for (i = 0; i < d->size; i ++) {

      delta = d->model[i] - d->mean[i];
      
      d->mean[i] += delta/(double)(d->counter);
      d->variance[i] += delta * (d->model[i] - d->mean[i]);
    }

    /*
     * Update the histogram
     */
    for (i = 0; i < d->size; i ++) {
      hi = histogram_index(d->model[i], d->vmin, d->vmax, d->bins);

      d->hist[i][hi] ++;
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

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int m;
  int mi;

  m = 0;
  mi = -1;

  for (i = 0; i < bins; i ++) {
    if (hist[i] > m) {
      m = hist[i];
      mi = i;
    }
  }
  
  if (mi < 0) {
    return 0.0;
  }

  return ((double)mi + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double median_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int j;
  int ci;
  int cj;

  i = 0;
  j = bins - 1;
  ci = 0;
  cj = 0;

  while (i != j) {
    if (ci < cj) {
      ci += hist[i];
      i ++;
    } else {
      cj += hist[j];
      j --;
    }
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = 0; 
  ci = 0;
  while(i < bins && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i ++;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = bins - 1; 
  ci = 0;
  while(i > 0 && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i --;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double hpd_from_histogram(int *hist, double vmin, double vmax, int bins, double hpd_interval, double &hpd_min, double &hpd_max)
{
  int sum;
  int mincount;
  
  //
  // First count number of samples
  //

  sum = 0;
  for (int i = 0; i < bins; i ++) {
    sum += hist[i];
  }

  mincount = (int)(hpd_interval * (double)sum);

  //
  // Now brute force search for minimum hpd with edges at
  //
  double minwidth = vmax - vmin;
  double minleft = vmin;
  double minright = vmax;
  
  for (int i = 0; i < bins; i ++) {

    double left = vmin + (double)i/(double)bins * (vmax - vmin);
    int j = i + 1;
    
    int count = hist[i];
    while (j < bins && count < mincount ) {
      count += hist[j];
      j ++;
    }

    if (count >= mincount) {

      double right = vmin + (double)j/(double)bins * (vmax - vmin);
      if (right - left < minwidth) {

	minwidth = right - left;
	minleft = left;
	minright = right;

      }

    }
  }

  hpd_min = minleft;
  hpd_max = minright;
  
  return minwidth;
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
	  " -o|--output <file>               Output mean model file\n"
	  " -v|--variance <file>             Output variance model file\n"
	  " -D|--stddev <file>               Output std dev. file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -m|--mode <file>                 Output mode\n"
	  " -M|--median <file>               Output median\n"
	  " -c|--credible-min <file>         Output credible min file\n"
	  " -C|--credible-max <file>         Output credible max file\n"
	  " -g|--histogram <file>            Output histogram file\n"
	  "\n"
	  " -b|--bins <int>                  No. histogram bins\n"
	  " -z|--vmin <float>                Lower range for histogram\n"
	  " -Z|--vmax <float>                Upper range for histogram\n"
	  "\n"
	  " -S|--maxsteps <int>              Chain history max steps\n"
	  "\n"
	  " -w|--wavelet-vertical <int>      Wavelet for vertical direction\n"
	  " -W|--wavelet-horizontal <int>    Wavelet for horizontal direction\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
