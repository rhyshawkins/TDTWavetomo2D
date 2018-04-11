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

#include "global.hpp"

extern "C" {
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "wavetree2d_sub.h"
}

static char short_options[] = "i:c:t:T:x:y:w:nZm:M:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"coefficients", required_argument, 0, 'c'},

  {"threshold", required_argument, 0, 't'},
  {"threshold-file", required_argument, 0, 'T'},
  
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},

  {"wavelet-xy", required_argument, 0, 'w'},

  {"norm", no_argument, 0, 'n'},

  {"zero-dc", no_argument, 0, 'Z'},
  {"mean", no_argument, 0, 'm'},
  {"mean-file", no_argument, 0, 'M'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void update_m4(double v, int *n, double *mean, double *meanabs, double *min, double *max);

static void usage(const char *pname);

static int load_image(const char *filename, int width, int height, double *img);
static int save_image(const char *filename, int width, int height, double *img);

constexpr int SUBTILE = 1;

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_model;
  char *coeff_file;

  double threshold;
  char *threshold_file;
  
  int degree_x;
  int degree_y;

  int width;
  int height;
  int workspace_size;
  int size;

  double *model;
  double *image;
  double *workspace;

  int i;
  int j;
  int k;
  int l;

  int degree_max;
  int ncoeff;
  
  double minc;
  double maxc;
  double meanc;
  double meanabsc;
  int meann;

  int wavelet_xy;

  generic_lift_forward1d_step_t waveletf;
  generic_lift_inverse1d_step_t waveleti;

  wavetree2d_sub_t *wt;

  bool norm;
  bool zerodc;

  double mean;

  char *mean_file;

  input_model = nullptr;
  coeff_file = nullptr;

  threshold = 0.1;
  threshold_file = nullptr;
  
  degree_x = 6;
  degree_y = 5;

  wavelet_xy = 0;

  norm = false;
  zerodc = false;

  mean = 0.0;
  mean_file = nullptr;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {

    case 'i':
      input_model = optarg;
      break;  

    case 'c':
      coeff_file = optarg;
      break;

    case 't':
      threshold = atof(optarg);
      if (threshold <= 0.0) {
	fprintf(stderr, "error: threshold must be greater than 0\n");
	return -1;
      }
      break;
      
    case 'T':
      threshold_file = optarg;
      break;

    case 'x':
      degree_x = atoi(optarg);
      if (degree_x < 1) {
	fprintf(stderr, "error: degree must be greater than 1\n");
	return -1;
      }
      break;

    case 'y':
      degree_y = atoi(optarg);
      if (degree_y < 1) {
	fprintf(stderr, "error: degree must be greater than 1\n");
	return -1;
      }
      break;
      
    case 'w':
      wavelet_xy = atoi(optarg);
      if (wavelet_xy < 0 || wavelet_xy > Global::WAVELET_MAX) {
	fprintf(stderr, "error: vertical wavelet must be between 0 and %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'n':
      norm = true;
      break;

    case 'Z':
      zerodc = true;
      break;

    case 'm':
      mean = atof(optarg);
      break;

    case 'M':
      mean_file = optarg;
      break;
      
    case 'h':
    default:
      usage(argv[0]);
      
      return -1;
    }
  }

  if (input_model == NULL) {
    fprintf(stderr, "error: require the input of a model file\n");
    return -1;
  }

  width = 1 << degree_x;
  height = 1 << degree_y;

  printf(" %d x %d degree\n", degree_x, degree_y);
  printf(" %d x %d image\n", width, height);
  
  size = width * height;
  model = new double[size];
  image = new double[size];

  if (load_image(input_model, width, height, model) < 0) {
    fprintf(stderr, "error: failed to load model\n");
    return -1;
  }

  workspace_size = width;
  if (height > workspace_size) {
    workspace_size = height;
  }
    
  workspace = new double[workspace_size];

  /*
   * Forward transform
   */
  waveletf = Global::wavelet_forward_function_from_id(wavelet_xy);
  waveleti = Global::wavelet_inverse_function_from_id(wavelet_xy);
  
  if (generic_lift_forward2d(model, 
			     width,
			     height,
			     width,
			     workspace,
			     waveletf,
			     waveletf,
			     SUBTILE) < 0) {
    fprintf(stderr, "error: failed to do forward transform\n");
    return -1;
  }

  if(save_image("test-coeff.txt", width, height, model) < 0) {
    fprintf(stderr, "error: faled to save full coeff\n");
    return -1;
  }

  for (int i = 0; i < size; i ++) {
    image[i] = model[i];
  }
  if (generic_lift_inverse2d(image,
			     width,
			     height,
			     width,
			     workspace,
			     waveleti,
			     waveleti,
			     SUBTILE) < 0) {
    fprintf(stderr, "error: failed to do inverse transform\n");
    return -1;
  }
  std::string fullimagefilename(threshold_file);
  fullimagefilename += ".fullimage";
  if (save_image(fullimagefilename.c_str(), width, height, image) < 0) {
    fprintf(stderr, "error: failed to save thresolded image\n");
    return -1;
  }
  
  wt = wavetree2d_sub_create(degree_x, degree_y, 0.0);
  if (wt == NULL) {
    fprintf(stderr, "error: failed to create wavetree\n");
    return -1;
  }

  printf(" %d x %d wavetree\n", wavetree2d_sub_get_width(wt), wavetree2d_sub_get_height(wt));
  
  ncoeff = wavetree2d_sub_get_ncoeff(wt);
  degree_max = wavetree2d_sub_maxdepth(wt);
  
  
  /*
   * Compute Coefficient statistics
   */
  minc = 1e9;
  maxc = -1e9;
  meanc = 0.0;
  meanabsc = 0.0;
  meann = 0;


  for (l = 1; l < ncoeff; l ++) {
    if (wavetree2d_sub_depthofindex(wt, l) == 1) {
      if (wavetree2d_sub_2dindices(wt, l, &i, &j) < 0) {
	fprintf(stderr, "error: failed to get 2d indices\n");
	return -1;
      }

      update_m4(model[j * width + i], &meann, &meanc, &meanabsc, &minc, &maxc);
    }
  }
  printf("%2d %10.6f %10.6f (%10.6f) %10.6f (%d)\n", 0, meanc, meanc, meanabsc, meanc, meann);
  
  for (k = 1; k <= degree_max; k ++) {

    minc = 1e9;
    maxc = -1e9;
    meanc = 0.0;
    meanabsc = 0.0;
    meann = 0;

    /* printf("%2d: %5d ... %5d\n", k + 1, degree_start, degree_end); */
    for (l = 1; l < ncoeff; l ++) {

      if (wavetree2d_sub_depthofindex(wt, l) == k) {
	c ++;
	if (wavetree2d_sub_2dindices(wt, l, &i, &j) < 0) {
	  fprintf(stderr, "error: failed to get 2d indices\n");
	  return -1;
	}
	
	update_m4(model[j * width + i], &meann, &meanc, &meanabsc, &minc, &maxc);
      }
    }
    
    printf("%2d %10.6f %10.6f (%10.6f) %10.6f (%d)\n", k, minc, meanc, meanabsc, maxc, meann);
  }

  if (coeff_file != NULL) {
    if (save_image(coeff_file, width, height, model) < 0) {
      fprintf(stderr, "error: failed to save coefficients\n");
      return -1;
    }
  }

  if (norm) {
    double l1norm = 0.0;
    for (int i = 1; i < size; i ++) {
      l1norm += fabs(model[i]);
    }
    printf("l1 %10.6f (%10.6f)\n", l1norm, l1norm/(double)(size - 1));
  }
	 
  if (threshold_file != nullptr) {
    if (mean_file != nullptr) {
      FILE *fp = fopen(mean_file, "r");
      if (fp == NULL) {
	fprintf(stderr, "error: failed to open mean file\n");
	return -1;
      }

      if (fscanf(fp, "%lf\n", &mean) != 1) {
	fprintf(stderr, "error: failed to read mean\n");
	return -1;
      }

      fclose(fp);
    }
    
    //
    // Create model
    //
    int base_size = wavetree2d_sub_base_size(wt);
    const int *base_indices = wavetree2d_sub_base_indices(wt);

    double modelmean = 0;

    if (base_size == 1) {
      modelmean = model[0];
    } else {
      for (int i = 0; i < base_size; i ++) {
	modelmean += model[base_indices[i] - 1];
      }
      modelmean /= (double)base_size;
    }
    printf("Model DC: %10.6f\n", modelmean);
    
    if (zerodc) {
      if (base_size == 1) {
	model[0] = 0.0;
      } else {
	for (int i = 0; i < base_size; i ++) {
	  model[base_indices[i] - 1] = 0.0;
	}
      }
      modelmean = 0.0;
    } else {
      if (base_size == 1) {
	model[0] -= mean;
	modelmean = model[0];
      } else {
	modelmean = 0.0;
	
	for (int i = 0; i < base_size; i ++) {
	  model[base_indices[i] - 1] -= mean;
	  modelmean += model[base_indices[i] - 1];
	}
	modelmean /= (double)base_size;
      }
    }

    printf("Model DC (demean): %10.6f\n", modelmean);
    
    if (wavetree2d_sub_create_from_array_with_threshold(wt, model, size, threshold) < 0) {
      fprintf(stderr, "error: failed to create thresholded model\n");
      return -1;
    }

    printf("Threshold: %.9g %d coeff\n", threshold, wavetree2d_sub_coeff_count(wt));

    //
    // Save wavetree model
    //
    if (wavetree2d_sub_save(wt, threshold_file) < 0) {
      fprintf(stderr, "error: failed to save thresholded model\n");
      return -1;
    }

    //
    // Reconstruct image.
    //
    memset(model, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, model, size) < 0) {
      fprintf(stderr, "error: failed to map thresholded model to array\n");
      return -1;
    }

    if(save_image("test-coeff-thresh.txt", width, height, model) < 0) {
      fprintf(stderr, "error: faled to save thresh coeff\n");
      return -1;
    }
    

    if (generic_lift_inverse2d(model, 
			       width,
			       height,
			       width,
			       workspace,
			       waveleti,
			       waveleti,
			       SUBTILE) < 0) {
      fprintf(stderr, "error: failed to do inverse transform\n");
      return -1;
    }

    if (!zerodc) {
      for (int i = 0; i < size; i ++) {
      	model[i] += mean;
      }
    }

    std::string filename(threshold_file);
    filename += ".image";
    if (save_image(filename.c_str(), width, height, model) < 0) {
      fprintf(stderr, "error: failed to save thresolded image\n");
      return -1;
    }
    
  }

  return 0;

}

static void update_m4(double v, int *n, double *mean, double *meanabs, double *min, double *max)
{
  double delta;
  
  (*n) ++;
  delta = v - (*mean);
  (*mean) += delta/(double)(*n);

  if (v < (*min)) {
    *min = v;
  }
  if (v > (*max)) {
    *max = v;
  }

  v = fabs(v);
  delta = v - (*meanabs);
  (*meanabs) += delta/(double)(*n);
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i | --input <filename>           Model raw image filename (required)\n"
	  " -c | --coefficients <filename>    Output raw coefficents (opt.)\n"
	  "\n"
	  " -t | --threshold <float>          Threshold value for thresholded model output\n"
	  " -T | --threshold-file <filename>  Threshold model output file (image written to filename.image)\n"
	  "\n"
	  " -x | --degree-x <int>             No. depth layers as power of 2\n"
	  " -y | --degree-y <int>             No. horizontal points as power of 2\n"
	  "\n"
	  " -w | --wavelet-xy <int>           Wavelet to use\n"
	  "\n"
	  " -Z | --zerodc                     Set the 0,0 component to zero\n"
	  " -m | --mean <float>               Subtract value from 0,0 component\n"
	  "\n"
	  " -n | --norm                       Print l1 norm of wavelet coefficients\n"
	  "\n"
	  " -h | --help              Show usage\n"
	  "\n",
	  pname);
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

static int
save_image(const char *filename,
           int width,
           int height,
           double *img)
{
  FILE *fp;
  int i;
  int j;
  double v;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "save_image: failed to create %s\n", filename);
    return -1;
  }

  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++) {
      v = img[j*width + i];
      fprintf(fp, "%.9g ", v);
    }
   fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}
