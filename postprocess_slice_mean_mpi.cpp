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

#include "globalslice.hpp"
#include "wavetomo2dutil.hpp"

struct user_data {
  int thincounter;
  int thin;
  int skip;
  
  int counter;
  
  int degree_max;
  int width;
  int height;
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

  wavetree2d_sub_t *wt;

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

static char short_options[] = "x:y:i:Z:o:D:t:s:m:M:c:C:g:p:P:Q:b:v:V:S:w:h";
static struct option long_options[] = {
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},

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

  int credible_drop;

  int waveleth;

  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
  /*
   * Default values
   */
  fp_in = NULL;
  fp_out = NULL;
  degree_x = 4;
  degree_y = 4;
  
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
  vmin = 1.0;
  vmax = 5.0;
  
  thin = 0;
  skip = 0;

  maxsteps = 1000000;

  waveleth = 0;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'x':
      degree_x = atoi(optarg);
      if (degree_x < 1) {
	fprintf(stderr, "error: invalid degree x\n");
	return -1;
      }
      break;

    case 'y':
      degree_y = atoi(optarg);
      if (degree_y < 1) {
	fprintf(stderr, "error: invalid degree y\n");
	return -1;
      }
      break;

    case 'i':
      input_file = optarg;
      break;

    case 'Z':
      zoffset_file = optarg;
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
      waveleth = atoi(optarg);
      if (waveleth < 0 || waveleth > GlobalSlice::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be between 0 and %d\n", (int)GlobalSlice::WAVELET_MAX);
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

  data.wt = wavetree2d_sub_create(degree_x, degree_y, 0.0);
  if (data.wt == nullptr) {
    fprintf(stderr, "error: failed to create wavetree\n");
    return -1;
  }
  
  data.width = wavetree2d_sub_get_width(data.wt);
  data.height = wavetree2d_sub_get_height(data.wt);
  data.size = wavetree2d_sub_get_size(data.wt);

  printf("Image: %d x %d\n", data.width, data.height);

  if (zoffset_file != NULL) {
    data.zoffset = load_zoffset(zoffset_file, data.zoffset_n);
    if (data.zoffset == NULL) {
      fprintf(stderr, "error: failed to load zoffset file\n");
      return -1;
    }

    if (data.zoffset_n != 1) {
      fprintf(stderr, "error: invalid zoffset file\n");
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
  for (int i = 0; i < data.size; i ++) {
    data.hist[i] = new int[data.bins];
    memset(data.hist[i], 0, sizeof(int) * data.bins);
  }
  
  data.variance = new double[data.size];
  memset(data.variance, 0, sizeof(double) * data.size);
  
  data.model = new double[data.size];

  data.workspace = new double[data.size];

  data.hwaveletf = GlobalSlice::wavelet_inverse_function_from_id(waveleth);

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }
  
  std::string chfile = mkfilenamerank(nullptr, input_file, mpi_rank);
  fp_in = fopen(chfile.c_str(), "r");
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

  for (int i = 0; i < data.size; i ++) {
    data.variance[i] /= (double)(data.counter - 1);
    data.variance[i] += (data.mean[i] * data.mean[i]);
  }

  /*
   * Mean output
   */
  MPI_Reduce(data.mean, data.workspace, data.size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    for (int i = 0; i < data.size; i ++) {
      data.mean[i] = data.workspace[i]/(double)mpi_size;
    }

    fp_out = fopen(output_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open mean file\n");
      return -1;
    }
    
    for (int j = 0; j < data.height; j ++) {
      for (int i = 0; i < data.width; i ++) {
	fprintf(fp_out, "%10.6f ", data.mean[j*data.width + i]);
      }
      fprintf(fp_out, "\n");
    }
    
    fclose(fp_out);
  }
 
  /*
   * Variance output
   */
  MPI_Reduce(data.variance, data.workspace, data.size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (mpi_rank == 0) {
    for (int i = 0; i < data.size; i ++) {
      data.variance[i] = data.workspace[i]/(double)mpi_size - data.mean[i]*data.mean[i];
    }

    /*
     * Std. Deviation output
     */
    if (stddev_file != NULL) {
      fp_out = fopen(stddev_file, "w");
      if (fp_out == NULL) {
	fprintf(stderr, "error: failed to open stddev file\n");
	return -1;
      }
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
	  fprintf(fp_out, "%10.6f ", sqrt(data.variance[j*data.width + i]));
	}
	fprintf(fp_out, "\n");
      }
      
      fclose(fp_out);
    }
  }

  if (mpi_rank == 0) {
    for (int i = 0; i < data.size; i ++) {
      MPI_Reduce(MPI_IN_PLACE, data.hist[i], data.bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  } else {
    for (int i = 0; i < data.size; i ++) {
      MPI_Reduce(data.hist[i], NULL, data.bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }

  if (mpi_rank == 0) {
    /*
     * Mode output
     */
    if (mode_file != NULL) {
      fp_out = fopen(mode_file, "w");
      if (fp_out == NULL) {
	fprintf(stderr, "error: failed to open mode file\n");
	return -1;
      }
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
	  fprintf(fp_out, "%10.6f ", mode_from_histogram(data.hist[j*data.width + i],
							 data.vmin, data.vmax, data.bins));
      }
	fprintf(fp_out, "\n");
      }
      
      fclose(fp_out);
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
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
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
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
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
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
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
      
      for (int j = 0; j < data.size; j ++) {
	for (int i = 0; i < data.bins; i ++) {
	  
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
      
      for (int j = 0; j < data.height; j ++) {
	for (int i = 0; i < data.width; i ++) {
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
  }
  
  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  delete [] data.mean;
  delete [] data.variance;
  delete [] data.model;
  delete [] data.workspace;

  MPI_Finalize();
  
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
  int hi;
  
  if ((d->thincounter >= d->skip) && (d->thin <= 1 || (d->thincounter % d->thin) == 0)) {
    
    memset(d->model, 0, sizeof(double) * d->size);

    if (wavetree2d_sub_set_from_S_v(d->wt, S_v) < 0) {
      fprintf(stderr, "process: failed to set wavetree (sub)\n");
      return -1;
    }

    if (wavetree2d_sub_map_to_array(d->wt, d->model, d->size) < 0) {
      fprintf(stderr, "process: failed to map to array\n");
      return -1;
    }

    
    if (generic_lift_inverse2d(d->model,
			       d->width,
			       d->height,
			       d->width,
			       d->workspace,
			       d->hwaveletf,
			       d->hwaveletf,
			       1) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    //
    // Add zoffset
    //
    if (d->zoffset != NULL) {
      for (int i = 0; i < d->size; i ++) {
	d->model[i] += d->zoffset[0];
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
	  " -x|--degree-x <int>              Number of longitude/x samples as power of 2\n"
	  " -y|--degree-y <int>              Number of latitude/y samples as power of 2\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  " -o|--output <file>               Output mean model file\n"
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
	  " -v|--vmin <float>                Lower range for histogram\n"
	  " -V|--vmax <float>                Upper range for histogram\n"
	  "\n"
	  " -S|--maxsteps <int>              Chain history max steps\n"
	  "\n"
	  " -w|--wavelet-lateral <int>       Wavelet for vertical direction\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
