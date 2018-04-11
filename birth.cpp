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

#include "birth.hpp"

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

Birth::Birth(Global &_global) :
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

Birth::~Birth()
{
  delete [] propose_depth;
  delete [] accept_depth;
}

int
Birth::step()
{
  propose ++;
  int k = wavetree3d_sub_coeff_count(global.wt);

  if (wavetree3d_sub_set_invalid_perturbation(global.wt, WT_PERTURB_BIRTH) < 0) {
    ERROR("failed to initialise birth perturbation\n");
    return -1;
  }

  if (k < global.kmax) {

    double ratio;
    int birth_depth;
    int birth_idx;
    double choose_prob;
    double birth_value;
    double birth_prob;
    int birth_valid = 0;
    double birth_parent_coeff;
    double reverse_prob;
    double prior_prob = 0.0;
    double proposed_likelihood;
    double proposed_log_normalization;
    int ii;
    int ij;
    int ik;
    
    if (choose_birth_location_and_value(k,
					ratio,
					birth_depth,
					birth_idx,
					choose_prob,
					birth_value,
					birth_prob,
					birth_valid,
					birth_parent_coeff,
					ii,
					ij,
					ik) < 0) {
      return -1;
    }

    if (communicate_birth_location_and_value(birth_valid,
					     birth_idx,
					     birth_depth,
					     birth_value) < 0) {
      return -1;
    }


    if (birth_valid) {

      propose_depth[birth_depth] ++;
      
      if (propose_birth(birth_valid,
			birth_idx,
			birth_depth,
			birth_value) < 0) {
	return -1;
      }
      
      if (compute_reverse_birth_probability(birth_depth,
					    birth_idx,
					    ii,
					    ij,
					    ik,
					    birth_parent_coeff,
					    birth_value,
					    reverse_prob,
					    prior_prob) < 0) {
	return -1;
      }
      
      if (compute_likelihood(birth_idx, proposed_likelihood, proposed_log_normalization) < 0) {
	return -1;
      }
    
      bool accept_proposal = false;

      if (compute_acceptance(proposed_likelihood,
			     proposed_log_normalization,
			     reverse_prob,
			     choose_prob,
			     birth_prob,
			     ratio,
			     prior_prob,
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
	accept_depth[birth_depth] ++;
	
	if (coefficient_histogram_accept_birth(global.coeff_hist, birth_idx, birth_value) < 0) {
	  ERROR("failed to update histogram for birth acceptance\n");
	  return -1;
	}

	if (wavetree3d_sub_commit(global.wt) < 0) {
	  ERROR("failed to commit birth\n");
	  return -1;
	}
	
	global.current_likelihood = proposed_likelihood;
	global.current_log_normalization = proposed_log_normalization;
	global.accept();
	
	return 1;
	
      } else {
	/*
	 * Reject
	 */
	
	if (coefficient_histogram_reject_birth(global.coeff_hist, birth_idx, birth_value) < 0) {
	  ERROR("failed to update histogram for birth rejection\n");
	  return -1;
	}
	
	if (wavetree3d_sub_undo(global.wt) < 0) {
	  ERROR("failed to undo birth\n");
	  return -1;
	}

	global.reject();
  
	return 0;
      }
    }
  }
  
  return 0;
}

std::string
Birth::write_short_stats()
{
  return mkformatstring("Birth %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
Birth::write_long_stats()
{
  std::string s = mkformatstring("Birth: %6d %7.3f:",
				 propose,
				 propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
  for (int i = 0; i <= global.treemaxdepth; i ++) {
    s = s + mkformatstring("%7.3f ",
			   propose_depth[i] == 0 ? 0.0 : 100.0*(double)accept_depth[i]/(double)propose_depth[i]);
  }
  
  return s;
}

void
Birth::initialize_mpi(MPI_Comm _communicator)
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
Birth::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
Birth::choose_birth_location_and_value(int k,
				       double &ratio,
				       int &birth_depth,
				       int &birth_idx,
				       double &choose_prob,
				       double &birth_value,
				       double &birth_prob,
				       int &birth_valid,
				       double &birth_parent_coeff,
				       int &ii,
				       int &ij,
				       int &ik)
{
  //
  // First determine coefficient to birth and new value
  //
  if (primary()) {

    if (hnk_get_kplus1_ratio(global.hnk, 
			     global.treemaxdepth, 
			     k, 
			     &ratio) < 0) {
      ERROR("failed to get ratio for birth\n");
      return -1;
    }
    
    if (ratio <= 0.0) {
      ERROR("invalid ratio %d %d %f\n", global.treemaxdepth, k, ratio);
      return -1;
    }
    
    //
    // Note for birth we actually want the inverse of the ratio
    //
    ratio = 1.0/ratio;
    
    if (wavetree3d_sub_choose_birth_global(global.wt,
					   global.random.uniform(),
					   global.treemaxdepth,
					   &birth_depth,
					   &birth_idx,
					   &choose_prob) < 0) {
      /* Generally means full tree, might need to check this later. */
      ERROR("failed to choose birth\n");
      return true;
    }
      
    if (coefficient_histogram_propose_birth(global.coeff_hist, birth_idx) < 0) {
      ERROR("failed to update histogram for birth proposal\n");
      return -1;
    }
      
    if (wavetree3d_sub_3dindices(global.wt, birth_idx, &ii, &ij, &ik) < 0) {
      ERROR("failed to compute 2d indices for birth\n");
      return -1;
    }
    
    if (wavetree3d_sub_get_coeff(global.wt,
				 wavetree3d_sub_parent_index(global.wt, birth_idx),
				 birth_depth - 1,
				 &birth_parent_coeff) < 0) {
      ERROR("failed to get parent coefficient for birth (idx = %d)\n", birth_idx);
      return -1;
    }
    
    if (wavetree_pp_birth3d(global.proposal,
			    ii, 
			    ij,
			    ik,
			    birth_depth,
			    global.treemaxdepth,
			    birth_parent_coeff,
			    &birth_value,
			    &birth_prob,
			    &birth_valid) < 0) {
      ERROR("Birth::step: failed to do birth proposal\n");
      return -1;
    }
  }

  return 0;
}

int
Birth::communicate_birth_location_and_value(int &birth_valid,
					    int &birth_idx,
					    int &birth_depth,
					    double &birth_value)
{
  if (communicator != MPI_COMM_NULL) {

    if (MPI_Bcast(&birth_valid, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast birth valid\n");
    }
    
    if (birth_valid) {
      if (MPI_Bcast(&birth_idx, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast birth index\n");
      }

      if (MPI_Bcast(&birth_depth, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast birth depth\n");
      }

      if (MPI_Bcast(&birth_value, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast birth value\n");
      }
    }
  }

  return 0;
}

int
Birth::propose_birth(int &birth_valid,
		     int &birth_idx,
		     int &birth_depth,
		     double &birth_value)
{
  if (birth_valid) {
    
    //
    // Birth the point (use birth from prior)
    //
    if (wavetree3d_sub_propose_birth(global.wt, 
				     birth_idx, 
				     birth_depth,
				     birth_value) < 0) {
      ERROR("failed to birth point\n");
      return -1;
    }

  }

  return 0;
}

int
Birth::compute_reverse_birth_probability(int birth_depth,
					 int birth_idx,
					 int ii,
					 int ij,
					 int ik,
					 double birth_parent_coeff,
					 double birth_value,
					 double &reverse_prob,
					 double &prior_prob)
{
  if (primary()) {
    if (wavetree3d_sub_reverse_birth_global(global.wt, 
					    global.treemaxdepth,
					    birth_depth,
					    birth_idx, 
					    &reverse_prob) < 0) {
      ERROR("failed to reverse birth (global)\n");
      return -1;
    }
    
    //
    // Compute the prior ratio
    //
    prior_prob = wavetree_pp_prior_probability3d(global.proposal,
						 ii, 
						 ij,
						 ik,
						 birth_depth,
						 global.treemaxdepth,
						 birth_parent_coeff,
						 birth_value);
  }

  return 0;
}

int
Birth::compute_likelihood(int birth_idx,
			  double &proposed_likelihood,
			  double &proposed_log_normalization)
{
  if (communicator == MPI_COMM_NULL) {
    proposed_likelihood = global.likelihood(proposed_log_normalization);
  } else {
    proposed_likelihood = global.likelihood_mpi(proposed_log_normalization);
  }

  return 0;
}

int
Birth::compute_acceptance(double proposed_likelihood,
			  double proposed_log_normalization,
			  double reverse_prob,
			  double choose_prob,
			  double birth_prob,
			  double ratio,
			  double prior_prob,
			  bool &accept_proposal)
{
  if (primary()) {
    double u = log(global.random.uniform());
    
    accept_proposal = u < ((global.current_likelihood - proposed_likelihood)/global.temperature /* Likelihood ratio */
			   + log(reverse_prob) - 
			   log(choose_prob)                                      /* Depth/Node proposal ratio */
			   - log(birth_prob)                                     /* Coefficient proposal */
			   + log(ratio)                                          /* Tree Prior */
			   + log(prior_prob)                                     /* Coefficient prior */
			   );

  }

  return 0;
}
  
int
Birth::communicate_acceptance(bool &accept_proposal)
{
  
  if (communicator != MPI_COMM_NULL) {

    int ta;

    if (mpi_rank == 0) {
      ta = (int)accept_proposal;
    }

    if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast accept proposal\n");
    }

    if (mpi_rank != 0) {
      accept_proposal = (bool)ta;
    }

  }

  return 0;
}

