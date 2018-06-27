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
#include <math.h>

#include <getopt.h>

#include <map>

#include "global.hpp"

#include "coordinate.hpp"
#include "wavetomo2dobservations.hpp"
#include "hierarchicalmodel.hpp"
#include "rng.hpp"

double synthetic_uniform3km(double nx, double ny, double nf)
{
  return 3.0;
}

double synthetic_cosinechecker3km_0(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI - M_PI) *
		      cos(ny * 2.0 * M_PI - M_PI) *
		      cos(nf * 2.0 * M_PI - M_PI));
}

double synthetic_cosinechecker3km_1(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 2.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 2.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 2.0 - M_PI));
}

double synthetic_cosinechecker3km_2(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 3.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 3.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 3.0 - M_PI));
}

double synthetic_cosinechecker3km_3(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 4.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 4.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 4.0 - M_PI));
}

double synthetic_cosinechecker3km_10(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 2.0 - M_PI) *
		      cos(ny * 2.0 * M_PI - M_PI) *
		      cos(nf * 2.0 * M_PI - M_PI));
}

double synthetic_cosinechecker3km_21(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 3.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 2.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 2.0 - M_PI));
}

double synthetic_cosinechecker3km_32(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 4.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 3.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 3.0 - M_PI));
}

double synthetic_cosinechecker3km_43(double nx, double ny, double nf)
{
  return 3.0 + 0.5 * (cos(nx * 2.0 * M_PI * 5.0 - M_PI) *
		      cos(ny * 2.0 * M_PI * 4.0 - M_PI) *
		      cos(nf * 2.0 * M_PI * 4.0 - M_PI));
}

double synthetic_ramplongitude3km(double nx, double ny, double nf)
{
  return 2.5 + nx;
}

double synthetic_ramplatitude3km(double nx, double ny, double nf)
{
  return 2.5 + ny;
}

double synthetic_rampfrequency3km(double nx, double ny, double nf)
{
  return 2.5 + nf;
}

static int image_width = -1;
static int image_height = -1;
static double *image = nullptr;

double synthetic_image(double nx, double ny, double nf)
{
  if (image != nullptr) {
    int i = nx * (double)image_width;
    int j = ny * (double)image_height;

    return image[j * image_width + i];
  }

  return -1.0;
}

typedef double (*synthetic_t)(double, double, double);
static std::map<std::string, synthetic_t> synthetic_models = {
  {"Uniform", synthetic_uniform3km},
  {"CosineCheck0", synthetic_cosinechecker3km_0},
  {"CosineCheck1", synthetic_cosinechecker3km_1},
  {"CosineCheck2", synthetic_cosinechecker3km_2},
  {"CosineCheck3", synthetic_cosinechecker3km_3},
  {"CosineCheck10", synthetic_cosinechecker3km_10},
  {"CosineCheck21", synthetic_cosinechecker3km_21},
  {"CosineCheck32", synthetic_cosinechecker3km_32},
  {"CosineCheck43", synthetic_cosinechecker3km_43},
  {"RampLongitude", synthetic_ramplongitude3km},
  {"RampLatitude", synthetic_ramplatitude3km},
  {"RampFrequency", synthetic_rampfrequency3km}
};

static char short_options[] = "o:t:I:U:T:m:M:x:y:z:u:n:N:a:A:f:F:Cs:W:H:S:lh";
static struct option long_options[] = {
  {"output", required_argument, 0, 'o'},
  {"template-input", required_argument, 0, 't'},
  {"image-output", required_argument, 0, 'I'},
  {"upscale-output", required_argument, 0, 'U'},
  {"true-output", required_argument, 0, 'T'},
  {"model", required_argument, 0, 'm'},
  {"model-image", required_argument, 0, 'M'},
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},
  {"super-resolution", required_argument, 0, 'u'},
  {"lonmin", required_argument, 0, 'n'},
  {"lonmax", required_argument, 0, 'N'},
  {"latmin", required_argument, 0, 'a'},
  {"latmax", required_argument, 0, 'A'},
  {"fmin", required_argument, 0, 'f'},
  {"fmax", required_argument, 0, 'F'},
  {"constant-z", no_argument, 0, 'C'},
  {"sigma", required_argument, 0, 's'},
  {"stations-x", required_argument, 0, 'W'},
  {"stations-y", required_argument, 0, 'H'},
  {"seed", required_argument, 0, 'S'},

  {"list-models", no_argument, 0, 'l'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static bool save_image(const char *filename, double *model, int width, int height);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Parameters
  //
  std::string modelname;
  const char *modelimage;
  
  int degreex;
  int degreey;
  int degreez;
  int super_resolution;

  double lonmin;
  double lonmax;
  double latmin;
  double latmax;

  double fmin;
  double fmax;

  double sigma;

  int nstationslon;
  int nstationslat;

  int seed;

  char *output;
  char *template_input;
  char *image_output;
  char *upscale_output;
  char *true_output;

  bool constant_z;

  bool list_models;

  int wavelet;

  //
  // Defaults
  //
  modelname = "Uniform";
  modelimage = nullptr;
  
  degreex = 4;
  degreey = 4;
  degreez = 4;
  super_resolution = 0;
  
  lonmin = -10.0;
  lonmax = 10.0;
  latmin = -10.0;
  latmax = 10.0;

  fmin = 0.3;
  fmax = 0.25;

  nstationslon = 4;
  nstationslat = 3;

  sigma = 0.1;

  seed = 983;

  output = nullptr;
  template_input = nullptr;
  image_output = nullptr;
  upscale_output = nullptr;
  true_output = nullptr;
  
  list_models = false;

  constant_z = false;

  wavelet = Global::WAVELET_CDF97;

  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'o':
      output = optarg;
      break;

    case 't':
      template_input = optarg;
      break;

    case 'I':
      image_output = optarg;
      break;

    case 'U':
      upscale_output = optarg;
      break;

    case 'T':
      true_output = optarg;
      break;
      
    case 'm':
      modelname = optarg;
      break;

    case 'M':
      modelimage = optarg;
      break;

    case 'x':
      degreex = atoi(optarg);
      if (degreex < 1 || degreex > 16) {
	fprintf(stderr, "error: degreex out of range\n");
	return -1;
      }
      break;

    case 'y':
      degreey = atoi(optarg);
      if (degreey < 1 || degreey > 16) {
	fprintf(stderr, "error: degreey out of range\n");
	return -1;
      }
      break;

    case 'z':
      degreez = atoi(optarg);
      if (degreez < 0 || degreez > 16) {
	fprintf(stderr, "error: degreez out of range\n");
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

    case 'f':
      fmin = atof(optarg);
      break;

    case 'F':
      fmax = atof(optarg);
      break;

    case 'C':
      constant_z = true;
      break;

    case 's':
      sigma = atof(optarg);
      if (sigma <= 0.0) {
	fprintf(stderr, "error: sigma must be greater than 0\n");
	return -1;
      }
      break;

    case 'W':
      nstationslon = atoi(optarg);
      if (nstationslon <= 0) {
	fprintf(stderr, "error: need at least one station in x direction\n");
	return -1;
      }
      break;

    case 'H':
      nstationslat = atoi(optarg);
      if (nstationslat <= 0) {
	fprintf(stderr, "error: need at least one station in y direction\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'l':
      list_models = true;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (list_models) {

    for (auto &m : synthetic_models) {
      printf("%s\n", m.first.c_str());
    }

    return 0;
  }
    

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }
  
  if ((nstationslon * nstationslat) <= 1) {
    fprintf(stderr, "error: not enough stations\n");
    return -1;
  }

  synthetic_t model_func = nullptr;

  if (modelimage != nullptr) {
    image_width = 1 << degreex;
    image_height = 1 << degreey;
    image = new double[image_width * image_height];
    FILE *fp = fopen(modelimage, "r");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open %s for reading\n", modelimage);
      return -1;
    }

    for (int j = 0; j < image_height; j ++) {
      for (int i = 0; i < image_width; i ++) {
	if (fscanf(fp, "%lf", &image[j * image_width + i]) != 1) {
	  fprintf(stderr, "error: failed to read image\n");
	  return -1;
	}
      }
    }

    fclose(fp);

    model_func = synthetic_image;
  } else {
    std::map<std::string, synthetic_t>::iterator it = synthetic_models.find(modelname);
    if (it == synthetic_models.end()) {
      fprintf(stderr, "error: invalid model name %s\n", modelname.c_str());
      return -1;
    }
    model_func = it->second;
  }

  generic_lift_inverse1d_step_t waveletxy = Global::wavelet_inverse_function_from_id(wavelet);
  
  //
  // Build the model of width x height x frequencies
  //
  int width = 1 << degreex;
  int height = 1 << degreey;
  int depth = 1 << degreez;
  int size = width * height * depth;

  printf("Model size: %d %d %d: %d\n", width, height, depth, size);
  
  double *model = new double[size];

  wavetomo2dobservations<LonLat<>> *obs;
  bool traveltime = false;

  if (template_input == nullptr) {
    obs = new wavetomo2dobservations<LonLat<>>(degreex, degreey, degreez,
					       super_resolution,
					       lonmin, lonmax, latmin, latmax,
					       traveltime,
					       false, // Linear
					       waveletxy);
  
    for (int k = 0; k < depth; k ++) {
      double nf = ((double)k + 0.5)/(double)depth;
      double f = fmin + (fmax - fmin) * nf;
      
      obs->add_frequency(f);
      
      for (int j = 0; j < height; j ++) {
	
	double nlat = ((double)j + 0.5)/(double)height;
	for (int i = 0; i < width; i ++) {
	  
	  double nlon = ((double)i + 0.5)/(double)width;
	  
	  //
	  // Model generated in normalized coordinates (i.e. 0 .. 1 in lon/lat/freq range)
	  //
	  if (constant_z) {
	    model[k * width * height + j * width + i] = model_func(nlon, nlat, 0.5);
	  } else {
	    model[k * width * height + j * width + i] = model_func(nlon, nlat, nf);
	  }
	  
	}
      }
    }
  

    //
    // Construct the Synthetic Stations
    //
    int idx = 0;
    for (int j = 0; j < nstationslat; j ++) {
      
      double slat = latmin + (latmax - latmin) * ((double)j + 0.5)/(double)nstationslat;
      
      for (int i = 0; i < nstationslon; i ++) {
	
	double slon = lonmin + (lonmax - lonmin) * ((double)i + 0.5)/(double)nstationslon;
	
	char code[15];
	sprintf(code, "XXX%03d", idx);
	idx ++;
	
	obs->add_station(code, LonLat<>(slon, slat));
      }
    }
    
    //
    // Construct the Traces between all station pairs
    //
    int nstations = nstationslat * nstationslon;
    
    for (int i = 0; i < (nstations - 1); i ++) {
      char codeA[15];
      sprintf(codeA, "XXX%03d", i);
      LonLat<> locationA;
      
      
      for (int j = i + 1; j < nstations; j ++) {
	char codeB[15];
	sprintf(codeB, "XXX%03d", j);
	
	obs->add_trace(codeA, codeB, true);
	
      }
    }
  
    obs->dump_traces();
  } else {

    obs = new wavetomo2dobservations<LonLat<>>(degreex, degreey, degreez,
					       super_resolution,
					       lonmin, lonmax, latmin, latmax,
					       template_input,
					       traveltime,
					       false, // Linear
					       waveletxy);

    if (depth != obs->get_frequency_count()) {
      fprintf(stderr, "error: frequency count does not match depth: %d != %d\n",
	      depth,
	      (int)obs->get_frequency_count());
      return -1;
    }
    
    double minfreq = obs->get_frequency(0);
    double maxfreq = obs->get_frequency(obs->get_frequency_count() - 1);

      
    for (int k = 0; k < depth; k ++) {
      double f = obs->get_frequency(k);
      double nf = (f - minfreq)/(maxfreq - minfreq);
      
      for (int j = 0; j < height; j ++) {
	
	double nlat = ((double)j + 0.5)/(double)height;
	for (int i = 0; i < width; i ++) {
	  
	  double nlon = ((double)i + 0.5)/(double)width;
	  
	  //
	  // Model generated in normalized coordinates (i.e. 0 .. 1 in lon/lat/freq range)
	  //
	  if (constant_z) {
	    model[k * width * height + j * width + i] = model_func(nlon, nlat, 0.5);
	  } else {
	    model[k * width * height + j * width + i] = model_func(nlon, nlat, nf);
	  }
	  
	  if (model[k * width * height + j * width + i] == 0.0) {
	    fprintf(stderr, "0 velocity: %d %f %f (%d, %d)\n", k, f, nf, i, j);
	    return -1;
	  }
	  
	}
      }
    }

  }
  
  //
  // Compute travel times and observed velocities
  //
  obs->compute_predictions(model, 0.0);
  obs->copy_predictions();

  if (upscale_output != nullptr) {
    if (!obs->save_super_image(upscale_output)) {
      fprintf(stderr, "error: failed to save upscaled image\n");
      return -1;
    }
  }

  //
  // Add noise
  //
  independentgaussianhierarchicalmodel hmodel;
  hmodel.setparameter(0, sigma);

  Rng rng(seed);

  if (true_output != nullptr) {
    if (!obs->save(true_output)) {
      fprintf(stderr, "error: failed to save true output\n");
      return -1;
    }

    char filename[1024];
    sprintf(filename, "%s.velocities", true_output);
    if (!obs->save_velocities(filename)) {
      fprintf(stderr, "error: failed to save true output velocities\n");
      return -1;
    }
      
  }

    
  obs->add_noise(rng, hmodel);

  if (!obs->save(output)) {
    fprintf(stderr, "error: failed to save output\n");
    return -1;
  }

  if (image_output != nullptr) {

    for (int d = 0; d < depth; d ++) {

      char filename[1024];
      sprintf(filename, "%s-%03d", image_output, d);
      if (!save_image(filename, model + d*width*height, width, height)) {
	fprintf(stderr, "error: failed to save image\n");
	return -1;
      }
    }
  }
	
    
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -o|--output <filename>              Output file to write (required)\n"
	  " -I|--image-output <filename>        Image output file prefix to write (optional)\n"
	  "\n"
	  " -m|--model <name>                   Synthetic model name\n"
	  "\n"
	  " -x|--degree-x <int>                 Lon samples as power of 2\n"
	  " -y|--degree-y <int>                 Lat samples as power of 2\n"
	  " -z|--degree-z <int>                 Freq samples as power of 2\n"
	  " -u|--super-resolution <int>         Upscale amount for ray tracing\n"
	  "\n"
	  " -n|--lonmin <float>                 Longitude min\n"
	  " -N|--lonmax <float>                 Longitude max\n"
	  " -a|--latmin <float>                 Latitude min\n"
	  " -A|--latmax <float>                 Latitude max\n"
	  " -f|--fmin <float>                   Frequency min\n"
	  " -F|--fmax <float>                   Frequency max\n"
	  "\n"
	  " -s|--sigma <float>                  Std deviation of noise\n"
	  "\n"
	  " -W|--stations-x <int>               No. stations along longitude\n"
	  " -H|--stations-y <int>               No. stations along latitude\n"
	  "\n"
	  " -S|--seed<int>                      Random number seed\n"
	  "\n"
	  " -l|--list-models                    List available synthetic models and exit\n"
	  "\n"
	  " -h|--help                           Show usage and exit\n"
	  "\n",
	  pname);
}

static bool save_image(const char *filename, double *model, int width, int height)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%15.9f ", model[j * width + i]);
      
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return true;
}

