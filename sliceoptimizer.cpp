#include <string.h>
#include <math.h>

#include <getopt.h>

#include <gmp.h>

extern "C" {
#include "wavetree2d_sub.h"
#include "wavetreepp.h"
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "slog.h"
};

#include "globalslice.hpp"
#include "birthslice.hpp"
#include "deathslice.hpp"
#include "valueslice.hpp"
#include "hierarchicalslice.hpp"
#include "hierarchicalpriorslice.hpp"

#include "wavetomo2dutil.hpp"

static char short_options[] = "i:I:M:o:x:y:z:u:s:n:N:a:A:t:k:w:Eh";


static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"prior-file", required_argument, 0, 'M'},
  {"output", required_argument, 0, 'o'},

  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},
  {"super-resolution", required_argument, 0, 'u'},
  
  {"slice", required_argument, 0, 's'},

  {"lonmin", required_argument, 0, 'n'},
  {"lonmax", required_argument, 0, 'N'},
  {"latmin", required_argument, 0, 'a'},
  {"latmax", required_argument, 0, 'A'},

  {"total", required_argument, 0, 't'},

  {"kmax", required_argument, 0, 'k'},

  {"wavelet-lateral", required_argument, 0, 'w'},

  {"linear", no_argument, 0, 'E'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Parameters
  //
  char *input_obs;
  char *initial_model;
  char *prior_file;
  char *output_prefix;

  int degreex;
  int degreey;
  int degreez;
  int super_resolution;
  
  int slice;
  
  double lonmin;
  double lonmax;
  double latmin;
  double latmax;

  int total;
  int seed;

  double lambda;
  int kmax;
  int wavelet_xy;

  int verbosity;
  bool linear;
  
  //
  // Defaults
  //
  input_obs = nullptr;
  initial_model = nullptr;
  prior_file = nullptr;
  output_prefix = nullptr;

  degreex = 7;
  degreey = 6;
  degreez = 5;
  super_resolution = 0;
  
  slice = 0;
  lonmin = -10.0;
  lonmax = 10.0;
  latmin = -10.0;
  latmax = 10.0;

  total = 10000;
  seed = 983;

  lambda = 1.0;
  kmax = 100;
  
  wavelet_xy = 0;

  verbosity = 1000;
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
      input_obs = optarg;
      break;

    case 'I':
      initial_model = optarg;
      break;

    case 'M':
      prior_file = optarg;
      break;

    case 'o':
      output_prefix = optarg;
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
      if (degreey < 1 || degreey > 16) {
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

    case 't':
      total = atoi(optarg);
      if (total < 1) {
        fprintf(stderr, "error: total must be greater than 0\n");
        return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'l':
      lambda = atof(optarg);
      if (lambda <= 0.0) {
        fprintf(stderr, "error: lambda std dev must be greater than 0\n");
        return -1;
      }
      break;
      
    case 'k':
      kmax = atoi(optarg);
      if (kmax < 1) {
        fprintf(stderr, "error: kmax must be greater than 0\n");
        return -1;
      }
      break;

    case 'w':
      wavelet_xy = atoi(optarg);
      if (wavelet_xy < 0 || wavelet_xy > GlobalSlice::WAVELET_MAX) {
        fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)GlobalSlice::WAVELET_MAX);
        return -1;
      }
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
                     initial_model,
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
                     seed,
                     kmax,
                     lambda,
                     false, // posterior ie Likelihood = constant
                     wavelet_xy,
                     linear);

  return 0;
}


  
static void usage(const char *pname)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "where options is one or more of:\n"
          "\n"
          " -i|--input <file>               Input observations file\n"
          " -I|--initial <file>             Starting model file\n"
          " -M|--prior <file>               Prior/Proposal file\n"
          " -o|--output <path>              Output prefix for output files\n"
          "\n"
          " -x|--degree-x <int>             Number of samples in x/lon direction as power of 2\n"
          " -y|--degree-y <int>             Number of samples in y/lat direction as power of 2\n"
          " -z|--degree-z <int>             Number of frequency samples as power of 2\n"
          " -s|--slice <int>                Index of slice to invert\n"
          "\n"
          " -n|--lonmin <float>             Longitude min\n"
          " -N|--lonmax <float>             Longitude max\n"
          " -a|--latmin <float>             Latitude min\n"
          " -A|--latmax <float>             Latitude max\n"
          "\n"
          " -t|--total <int>                Total number of iterations\n"
          "\n"
          " -k|--kmax <int>                 Max. no. of coefficients\n"
          "\n"
          " -w|--wavelet-xy <int>           Wavelet basis to use for lon/lat plane\n"
          "\n"
          " -h|--help                       Show usage information\n"
          "\n",
          pname);
}

