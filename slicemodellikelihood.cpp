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

#include <getopt.h>

#include "globalslice.hpp"

extern "C" {
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "wavetree2d_sub.h"
}

//
// Todo: Threshold wavelet model from image + likelihood
//

static char short_options[] = "i:I:m:o:x:y:z:u:s:n:N:a:A:w:Z:W:l:r:p:Eh";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"image", required_argument, 0, 'I'},
  {"mean", required_argument, 0, 'm'},
  {"observations", required_argument, 0, 'o'},
  
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},
  {"super-resolution", required_argument, 0, 'u'},
  
  {"slice", required_argument, 0, 's'},
  
  {"lonmin", required_argument, 0, 'n'},
  {"lonmax", required_argument, 0, 'N'},
  {"latmin", required_argument, 0, 'a'},
  {"latmax", required_argument, 0, 'A'},

  {"wavelet-xy", required_argument, 0, 'w'},

  {"zoffset", required_argument, 0, 'Z'},
  {"write-zoffset", required_argument, 0, 'W'},

  {"lambda", required_argument, 0, 'l'},
  {"residual", required_argument, 0, 'r'},
  {"predictions", required_argument, 0, 'p'},
  {"linear", no_argument, 0, 'E'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static int load_image(const char *filename, int width, int height, double *img);
static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_model;
  char *input_image;
  double mean;
  char *input_obs;
  char *residuals_file;
  char *predictions_file;
  
  int degreex;
  int degreey;
  int degreez;
  int super_resolution;
  
  int slice;
  
  double lonmin;
  double lonmax;
  double latmin;
  double latmax;

  int wavelet_xy;

  double lambda_scale;

  double zoffset;
  char *write_zoffset;

  bool linear;

  //
  // Defaults
  //

  input_model = nullptr;
  input_image = nullptr;
  mean = -1.0;
  input_obs = nullptr;
  residuals_file = nullptr;
  predictions_file = nullptr;
  
  degreex = 6;
  degreey = 5;
  degreez = 5;
  super_resolution = 0;
  
  slice = 0;

  lonmin = -10.0;
  lonmax = 10.0;
  latmin = -10.0;
  latmax = 10.0;

  wavelet_xy = 0;

  lambda_scale = 1.0;
  zoffset = 0.0;
  write_zoffset = nullptr;

  linear = false;
  
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
      input_model = optarg;
      break;

    case 'I':
      input_image = optarg;
      break;

    case 'm':
      mean = atoi(optarg);
      if (mean <= 0.0) {
	fprintf(stderr, "error: mean must be greater than 0\n");
	return -1;
      }
      break;
      
    case 'o':
      input_obs = optarg;
      break;

    case 'x':
      degreex = atoi(optarg);
      if (degreex < 1 || degreex > 16) {
	fprintf(stderr, "error: degree x must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'y':
      degreey = atoi(optarg);
      if (degreey < 1 || degreey > 16) {
	fprintf(stderr, "error: degree y must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'z':
      degreez = atoi(optarg);
      if (degreez < 1 || degreez > 16) {
	fprintf(stderr, "error: degree z must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'u':
      super_resolution = atoi(optarg);
      if (super_resolution < 0) {
	fprintf(stderr, "error: super resolution must be 0 or greater\n");
	return -1;
      }
      break;

    case 's':
      slice = atoi(optarg);
      if (slice < 0) {
	fprintf(stderr, "error: slice must be 0 or greater\n");
	return -1;
      }
      break;

    case 'n':
      lonmin = atof(optarg);
      break;

    case 'N':
      lonmax = atof(optarg);
      break;

    case 'a':
      latmin = atof(optarg);
      break;

    case 'A':
      latmax = atof(optarg);
      break;

    case 'l':
      lambda_scale = atof(optarg);
      break;
      
    case 'w':
      wavelet_xy = atoi(optarg);
      if (wavelet_xy < 0 || wavelet_xy > GlobalSlice::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)GlobalSlice::WAVELET_MAX);
	return -1;
      }
      break;

    case 'Z':
      zoffset = atof(optarg);
      break;

    case 'W':
      write_zoffset = optarg;
      break;

    case 'r':
      residuals_file = optarg;
      break;

    case 'p':
      predictions_file = optarg;
      break;

    case 'E':
      linear = true;
      break;
      
    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (input_obs == nullptr) {
    fprintf(stderr, "error: required input parameter input observations missing\n");
    return -1;
  }

  GlobalSlice global(input_obs,
		     input_model,
		     nullptr,
		     degreex,
		     degreey,
		     degreez,
		     super_resolution,
		     slice,
		     lonmin,
		     lonmax,
		     latmin,
		     latmax,
		     0,
		     100,
		     lambda_scale,
		     false,
		     wavelet_xy,
		     linear);

  double *image_model = nullptr;
  if (input_image != nullptr || mean > 0.0) {
    int width = 1 << degreex;
    int height = 1 << degreey;
    int size = width * height;
    image_model = new double[size];

    if (input_image != nullptr) {
      if (load_image(input_image, width, height, image_model) < 0) {
	fprintf(stderr, "error: failed to load image\n");
	return -1;
      }

      if (zoffset != 0.0) {
	printf("Configured Zoffset: %g\n", zoffset);
	for (int i = 0; i < size; i ++) {
	  image_model[i] += zoffset;
	}
      }
	
    } else {
      for (int i = 0; i < size; i ++) {
	image_model[i] = mean;
      }
    }

    double log_normalization;
    double like = global.image_likelihood(image_model, log_normalization);
    printf("Likelihood: %g (%g)\n", like, log_normalization);
  } else {
    
    double log_normalization;
    double like = global.likelihood(log_normalization);
    printf("Zoffset: %g\n", global.zoffset[slice]);
    printf("Likelihood: %g (%g)\n", like, log_normalization);

    if (write_zoffset != nullptr) {
      FILE *fp = fopen(write_zoffset, "w");
      if (fp == NULL) {
	fprintf(stderr, "error: failed to create zoffset output file\n");
	return -1;
      }

      fprintf(fp, "%15.9f\n", global.zoffset[slice]);
      fclose(fp);
    }
  }

  global.accept();

  if (residuals_file) {
    int n = global.get_residual_size();
    const double *mean_residuals = global.get_mean_residuals();
    const double *mean_normed_residuals = global.get_mean_normed_residuals();

    if (!global.observations->save_slice_residuals(residuals_file,
						   n,
						   mean_residuals,
						   mean_normed_residuals)) {
      fprintf(stderr, "error: failed to create residual file\n");
      return -1;
    }
  }

  if (predictions_file) {
    if (!global.observations->save_slice_predictions(predictions_file, slice)) {
      fprintf(stderr, "error: failed to create predictions file\n");
      return -1;
    }
  }
      

  return 0;
}
  
static int load_image(const char *filename, int width, int height, double *img)
{
  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "load_image: failed to open file\n");
    return -1;
  }

  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++) {

      if (fscanf(fp, "%lf", &(img[j*width + i])) != 1) {
        fprintf(stderr, "load_image: failed to read pixel\n");
        return -1;
      }

    }
  }

  fclose(fp);
  return 0;
}

static void
usage(const char *pname)
{
  fprintf(stderr, "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  "-i|--input <filename>          Input wavetree model file\n"
	  "-I|--image <filename>          Input image file\n"
	  "-o|--observations <filename>   Input observations file\n"
	  "\n"
	  "-x|--degree-x <int>            No. layers as a power of 2\n"
	  "-y|--degree-y <int>            No. lateral points as power of 2\n"
	  "-z|--degree-z <int>            No. frequencies as power of 2\n"
	  "-s|--slice <int>               Index of slice to computing likelihood for\n"
	  "\n"
	  " -n|--lonmin <float>           Longitude min\n"
	  " -N|--lonmax <float>           Longitude max\n"
	  " -a|--latmin <float>           Latitude min\n"
	  " -A|--latmax <float>           Latitude max\n"
	  "\n"
	  "-w|--wavelet-xy <int>          Wavelet in xy-plane\n"
	  "\n"
	  "-l|--lambda-scale <float>      Hierarchical lambda scaling parameter\n"
	  "\n"
	  "-h|--help                      Usage\n"
	  "\n",
	  pname);
}
