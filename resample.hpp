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

#pragma once
#ifndef resample_hpp

#include <string>

#include "global.hpp"
#include <mpi.h>

class Resample {
public:

  Resample(Global &global);
  ~Resample();

  int step(double resample_temperature = 1.0);

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm global_communicator,
		      MPI_Comm temperature_communicator,
		      MPI_Comm chain_communicator);

  Global &global;

  MPI_Comm global_communicator;
  MPI_Comm temperature_communicator;
  MPI_Comm chain_communicator;

  int resamplings;
  int reselected;
  int propagated;
  int replaced;

  int global_size;
  int global_rank;

  int temperature_size;
  int temperature_rank;

  int chain_size;
  int chain_rank;

  int processesperchain;

  double *likelihood_array;
  double *lambda_array;
  double *weights_array;
  int *sources_array;
  MPI_Request *requests_array;

  int send_buffer_size;
  char *send_buffer;

  int recv_buffer_size;
  char *recv_buffer;

};

#endif // resample_hpp
