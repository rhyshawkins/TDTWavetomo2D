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
#ifndef hierarchicalslice_hpp
#define hierarchicalslice_hpp

#include <mpi.h>

#include "globalslice.hpp"

extern "C" {
  #include "chain_history.h"
};

class HierarchicalSlice {
public:

  HierarchicalSlice(GlobalSlice &global, double sigma);
  ~HierarchicalSlice();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  void get_last_step(chain_history_change_t *last_step);

  GlobalSlice &global;
  double sigma;

  int propose;
  int accept;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

  chain_history_change_t last_step;

private:

  bool primary() const;

  int choose_value(double &value,
		   double &value_prior_ratio,
		   int &valid_proposal);
  
  int communicate_value(int &valid_proposal,
			double &value);

  int compute_likelihood(double value,
			 double &proposed_likelihood,
			 double &proposed_log_normalization);

  int compute_acceptance(double log_value_prior_ratio,
			 double proposed_likelihood,
			 double proposed_log_normalization,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);

};

#endif // hierarchicalslice_hpp
