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
#ifndef value_hpp
#define value_hpp

#include <mpi.h>

#include "global.hpp"

class Value {
public:

  Value(Global &global);
  ~Value();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  Global &global;

  int propose;
  int accept;

  int *propose_depth;
  int *accept_depth;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

private:

  bool primary() const;

  int choose_value_location_and_value(int &value_depth,
				      int &value_idx,
				      double &choose_prob,
				      double &value,
				      int &ii,
				      int &ij,
				      int &ik,
				      double &value_prior_ratio,
				      int &prior_errors,
				      int &valid_proposal);

  int communicate_value_location_and_value(int &valid_proposal,
					   int &value_idx,
					   int &value_depth,
					   double &value);

  int propose_value(int valid_proposal,
		    int value_idx,
		    int value_depth,
		    double value);

  int compute_likelihood(int valid_idx, double &proposed_likelihood, double &proposed_log_normalization);

  int compute_acceptance(int value_idx,
			 double value_prior_ratio,
			 double proposed_likelihood,
			 double proposed_log_normalization,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);
};

#endif // value_hpp
