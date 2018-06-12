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

static char short_options[] = "i:o:x:y:z:w:n:N:a:A:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},

  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},

  {"lonmin", required_argument, 0, 'n'},
  {"lonmax", required_argument, 0, 'N'},
  {"latmin", required_argument, 0, 'a'},
  {"latmax", required_argument, 0, 'A'},

  {"wavelet", required_argument, 0, 'w'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static double *generate_wavelet_matrix(int degreex, int degreey, generic_lift_inverse1d_step_t basis);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Parameters
  //
  char *input_obs;
  char *output_file;
  
  int degreex;
  int degreey;
  int degreez;
  int slice;
  
  double lonmin;
  double lonmax;
  double latmin;
  double latmax;

  int wavelet;
  
  //
  // Defaults
  //

  input_obs = nullptr;
  output_file = nullptr;
  
  degreex = 7;
  degreey = 6;
  degreez = 5;
  slice = 0;

  lonmin = -10.0;
  lonmax = 10.0;
  latmin = -10.0;
  latmax = 10.0;

  wavelet = 0;
  
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

    case 'o':
      output_file = optarg;
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

    case 'w':
      wavelet = atoi(optarg);
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

  if (output_file == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  wavetomo2dobservations<LonLat<>> *observations = new 
    wavetomo2dobservations<LonLat<>>(degreex,
				     degreey,
				     degreez,
				     0,
				     lonmin,
				     lonmax,
				     latmin,
				     latmax,
				     input_obs,
				     false,
				     true,
				     cdf97_lift_inverse1d_cdf97_step);

  printf("Cfg  Longitude: %10.6f %10.6f\n", observations->get_xmin(), observations->get_xmax());
  printf("Cfg  Latitude : %10.6f %10.6f\n", observations->get_ymin(), observations->get_ymax());

  printf("Data Longitude: %10.6f %10.6f\n", observations->get_station_xmin(), observations->get_station_xmax());
  printf("Data Latitude : %10.6f %10.6f\n", observations->get_station_ymin(), observations->get_station_ymax());

  double RADIUS = wavetomo2dobservations<LonLat<>>::RADIUS;
  double hdisttop = LinearWeights::gc_dist(lonmin, latmax, lonmax, latmax, RADIUS);
  double hdistbot = LinearWeights::gc_dist(lonmin, latmin, lonmax, latmin, RADIUS);
  printf("Horizontal : %10.6f %10.6f\n", hdisttop, hdistbot);
  printf("Vertical   : %10.6f\n", LinearWeights::gc_dist(lonmin, latmin, lonmin, latmax, RADIUS));
  
  if (!observations->compute_linear_weights()) {
    fprintf(stderr, "error: failed to compute linear weights\n");
    return -1;
  }

  if (!observations->save_identity_projection_matrix(output_file)) {
    fprintf(stderr, "error: failed to save projetion matrix\n");
    return -1;
  }

  char filename[1024];
  sprintf(filename, "%s.b_dist_error", output_file);
  if (!observations->save_b_dist_error_matrix(filename, slice)) {
    fprintf(stderr, "error: failed to b/dist/error matrix\n");
    return -1;
  }

  delete observations;

  return 0;
};
  
static void usage(const char *pname)
{
}

static double *generate_wavelet_matrix(int degreex, int degreey, generic_lift_inverse1d_step_t basis)
{
  int imagewidth = 2 << degreex;
  int imageheight = 2 << degreey;
  
  int N = imagewidth * imageheight;
  int size = (N*N);

  double *matrix = new double[size];

  
}
