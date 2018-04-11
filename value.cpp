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

#include <math.h>

extern "C" {
#include "slog.h"
};

#include "value.hpp"

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

Value::Value(Global &_global) :
  global(_global),
  propose(0),
  accept(0),
  propose_depth(new int[global.treemaxdepth + 1]),
  accept_depth(new int[global.treemaxdepth + 1]),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1)
{
  for (int i = 0; i <= global.treemaxdepth; i ++) {
    propose_depth[i] = 0;
    accept_depth[i] = 0;
  }
}

Value::~Value()
{
  delete [] propose_depth;
  delete [] accept_depth;
}

int
Value::step()
{
  propose ++;

  if (wavetree3d_sub_set_invalid_perturbation(global.wt, WT_PERTURB_VALUE) < 0) {
    return -1;
  }

  int value_depth;
  int value_idx;
  double choose_prob;
  double value;
  int ii, ij, ik;
  double proposed_likelihood;
  double proposed_log_normalization;
  double value_prior_ratio = 1.0; // This must be 1 before calling
  int valid_proposal = 0;
  int prior_errors = 0;

  if (choose_value_location_and_value(value_depth,
				      value_idx,
				      choose_prob,
				      value,
				      ii,
				      ij,
				      ik,
				      value_prior_ratio,
				      prior_errors,
				      valid_proposal) < 0) {
    return -1;
  }
  
  if (communicate_value_location_and_value(valid_proposal,
					   value_idx,
					   value_depth,
					   value) < 0) {
    return -1;
  }
  
  if (valid_proposal) {

    propose_depth[value_depth] ++;

    if (propose_value(valid_proposal,
		      value_idx,
		      value_depth,
		      value) < 0) {
      return -1;
    }

    if (compute_likelihood(value_idx, proposed_likelihood, proposed_log_normalization) < 0) {
      return -1;
    }

    bool accept_proposal = false;
    if (compute_acceptance(value_idx,
			   value_prior_ratio,
			   proposed_likelihood,
			   proposed_log_normalization,
			   accept_proposal) < 0) {
      return -1;
    }

    if (communicate_acceptance(accept_proposal) < 0) {
      return -1;
    }
      
    if (accept_proposal) {

      //
      // Accept
      //

      accept ++;
      accept_depth[value_depth] ++;
      
      if (coefficient_histogram_accept_value(global.coeff_hist, value_idx, value) < 0) {
        ERROR("failed to update histogram for value acceptance\n");
        return -1;
      }

      if (wavetree3d_sub_commit(global.wt) < 0) {
        ERROR("failed to commit value (%d %d)\n", valid_proposal, prior_errors);
        return -1;
      }
      
      global.current_likelihood = proposed_likelihood;
      global.current_log_normalization = proposed_log_normalization;
      global.accept();
      
      return 1;
      
    } else {
      
      //
      // Reject
      //
      
      if (coefficient_histogram_reject_value(global.coeff_hist, value_idx, value) < 0) {
        ERROR("failed to update histogram for value rejection\n");
        return -1;
      }

      if (wavetree3d_sub_undo(global.wt) < 0) {
        ERROR("failed to undo value\n");
        return -1;
      }

      global.reject();

      return 0;
    }

  }
  
  return 0;
}

std::string
Value::write_short_stats()
{
  return mkformatstring("Value %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
Value::write_long_stats()
{
  std::string s = mkformatstring("Value: %6d %7.3f:",
			       propose,
			       propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);

  for (int i = 0; i <= global.treemaxdepth; i ++) {
    s = s + mkformatstring("%7.3f ",
			   propose_depth[i] == 0 ? 0.0 : 100.0*(double)accept_depth[i]/(double)propose_depth[i]);
  }

  return s;
}

void
Value::initialize_mpi(MPI_Comm _communicator)
{
  MPI_Comm_dup(_communicator, &communicator);

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }
}

bool
Value::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
Value::choose_value_location_and_value(int &value_depth,
				       int &value_idx,
				       double &choose_prob,
				       double &value,
				       int &ii,
				       int &ij,
				       int &ik,
				       double &value_prior_ratio,
				       int &prior_errors,
				       int &valid_proposal)
{
  if (primary()) {

    if (wavetree3d_sub_choose_value_global(global.wt,
					   global.random.uniform(),
					   global.treemaxdepth,
					   &value_depth,
					   &value_idx,
					   &choose_prob) < 0) {
      ERROR("failed to choose global value\n");
      return -1;
    }
    
    if (wavetree3d_sub_get_coeff(global.wt,
				 value_idx,
				 value_depth,
				 &value) < 0) {
      ERROR("failed to get coefficient value\n");
      return -1;
    }
    
    if (wavetree3d_sub_3dindices(global.wt, value_idx, &ii, &ij, &ik) < 0) {
      ERROR("failed to compute 2d indices for birth\n");
      return -1;
    }
    
    if (coefficient_histogram_propose_value(global.coeff_hist, value_idx) < 0) {
      ERROR("failed to update histogram for value proposal\n");
      return -1;
    }
    
    if (wavetree_pp_value_init(global.proposal) < 0) {
      ERROR("failed to initialize value proposal\n");
      return -1;
    }

    double value_parent_coeff = 0.0;
    
    if (wavetree_pp_propose_value3d(global.proposal, 
				    ii, ij, ik,
				    value_depth, 
				    global.maxdepth, 
				    value_parent_coeff,
				    global.temperature,
				    &value,
				    &value_prior_ratio) < 0) {
      ERROR("failed to perturb value\n");
      return -1;
    }
  
    prior_errors = wavetree_pp_value_error_count(global.proposal);
    if (prior_errors < 0) {
      ERROR("failed to check errors\n");
      return -1;
    }
    
    if (prior_errors == 0) {
      valid_proposal = 1;
    }
  }

  return 0;
}
  

int
Value::communicate_value_location_and_value(int &valid_proposal,
					    int &value_idx,
					    int &value_depth,
					    double &value)
{
  if (communicator != MPI_COMM_NULL) {

    if (MPI_Bcast(&valid_proposal, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast valid proposal\n");
    }
    
    if (valid_proposal) {

      if (MPI_Bcast(&value_idx, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast index\n");
      }
      if (MPI_Bcast(&value_depth, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast depth\n");
      }
      if (MPI_Bcast(&value, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast value\n");
      }
      
    }
  }

  return 0;
}
  
int
Value::propose_value(int valid_proposal,
		     int value_idx,
		     int value_depth,
		     double value)
{
  if (valid_proposal) {
    if (wavetree3d_sub_propose_value(global.wt, value_idx, value_depth, value) < 0) {
      ERROR("failed to propose value\n");
      return -1;
    }
  }

  return 0;
}

int
Value::compute_likelihood(int value_idx, double &proposed_likelihood, double &proposed_log_normalization)
{
  if (communicator == MPI_COMM_NULL) {
    proposed_likelihood = global.likelihood(proposed_log_normalization);
  } else {
    proposed_likelihood = global.likelihood_mpi(proposed_log_normalization);
  }

  return 0;
}

int
Value::compute_acceptance(int value_idx,
			  double value_prior_ratio,
			  double proposed_likelihood,
			  double proposed_log_normalization,
			  bool &accept_proposal)
{
  if (primary()) {
    
    double u = log(global.random.uniform());
    
    double alpha = (log(value_prior_ratio) + (global.current_likelihood - proposed_likelihood)/global.temperature);
    
    if (coefficient_histogram_sample_value_alpha(global.coeff_hist, value_idx, exp(alpha)) < 0) {
      ERROR("failed to sample alpha\n");
      return -1;
    }
    
    accept_proposal = (u < alpha);

  }

  return 0;
}

int
Value::communicate_acceptance(bool &accept_proposal)
{
  
  if (communicator != MPI_COMM_NULL) {

    int ta;

    if (mpi_rank == 0) {
      ta = (int)accept_proposal;
    }

    if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast acceptted\n");
    }

    if (mpi_rank != 0) {
      accept_proposal = (bool)ta;
    }

  }

  return 0;
}
