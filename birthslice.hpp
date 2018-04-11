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
#ifndef birthslice_hpp
#define birthslice_hpp

#include "globalslice.hpp"

#include <mpi.h>

class BirthSlice {
public:

  BirthSlice(GlobalSlice &global);
  ~BirthSlice();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  GlobalSlice &global;
  int propose;
  int accept;

  int *propose_depth;
  int *accept_depth;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

private:

  bool primary() const;

  int choose_birth_location_and_value(int k,
				      double &ratio,
				      int &birth_depth,
				      int &birth_idx,
				      double &choose_prob,
				      double &birth_value,
				      double &birth_prob,
				      int &birth_valid,
				      double &birth_parent_coeff,
				      int &ii,
				      int &ij);

  int communicate_birth_location_and_value(int &birth_valid,
					   int &birth_idx,
					   int &birth_depth,
					   double &birth_value);

  int propose_birth(int &birth_valid,
		    int &birth_idx,
		    int &birth_depth,
		    double &birth_value);

  int compute_reverse_birth_probability(int birth_depth,
					int birth_idx,
					int ii,
					int ij,
					double birth_parent_coeff,
					double birth_value,
					double &reverse_prob,
					double &prior_prob);

  int compute_likelihood(int birth_idx, double &proposed_likelihood, double &proposed_log_normalization);

  int compute_acceptance(double proposed_likelihood,
			 double proposed_log_normalization,
			 double reverse_prob,
			 double choose_prob,
			 double birth_prob,
			 double ratio,
			 double prior_prob,
			 bool &accept_proposal);

  int communicate_acceptance(bool &accept_proposal);

};

#endif // birthslice_hpp
