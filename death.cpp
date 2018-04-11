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

#include "death.hpp"

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

Death::Death(Global &_global) :
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

Death::~Death()
{
  delete [] propose_depth;
  delete [] accept_depth;
}

int
Death::step()
{
  propose ++;

  int k = wavetree3d_sub_coeff_count(global.wt);

  if (wavetree3d_sub_set_invalid_perturbation(global.wt, WT_PERTURB_DEATH) < 0) {
    return -1;
  }

  if (k > 1) {

    double ratio;
    int death_valid = 1;
    int death_depth;
    int death_idx;
    double choose_prob;
    double death_value;
    int ii, ij, ik;
    double death_parent_coeff;
    double death_prob;
    double reverse_prob;
    double prior_prob;
    double proposed_likelihood;
    double proposed_log_normalization;
    
    //
    // First determine coefficient to death
    //
    if (choose_death_location(k,
			      ratio,
			      death_depth,
			      death_idx,
			      choose_prob,
			      death_valid) < 0) {
      return -1;
    }

    if (communicate_death_location(death_valid,
				   death_idx,
				   death_depth) < 0) {
      return -1;
    }

    
    if (death_valid) {

      propose_depth[death_depth] ++;
      
      if (propose_death(death_valid,
			death_idx,
			death_depth,
			death_value) < 0) {
	return -1;
      }
      
      if (compute_reverse_death_probability(death_idx,
					    death_depth,
					    death_value,
					    ii,
					    ij,
					    ik,
					    death_parent_coeff,
					    death_prob,
					    reverse_prob,
					    prior_prob) < 0) {
	return -1;
      }

      if (compute_likelihood(death_idx,
			     proposed_likelihood,
			     proposed_log_normalization) < 0) {
	return -1;
      }

      bool accept_proposal = false;

      if (compute_acceptance(proposed_likelihood,
			     proposed_log_normalization,
			     reverse_prob,
			     choose_prob,
			     death_prob,
			     ratio,
			     prior_prob,
			     accept_proposal) < 0) {
	return -1;
      }

      if (communicate_acceptance(accept_proposal) < 0) {
	return -1;
      }
      
      /*
       * Accept
       */
      if (accept_proposal) {
	accept ++;
	accept_depth[death_depth] ++;
	
        if (coefficient_histogram_accept_death(global.coeff_hist, death_idx) < 0) {
          ERROR("failed to update histogram for death acceptance\n");
          return -1;
        }
	
        if (wavetree3d_sub_commit(global.wt) < 0) {
          ERROR("failed to commit death\n");
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
        if (wavetree3d_sub_undo(global.wt) < 0) {
          ERROR("failed to undo death\n");
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
Death::write_short_stats()
{
  return mkformatstring("Death %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
Death::write_long_stats()
{
  std::string s = mkformatstring("Death: %6d %7.3f:",
				 propose,
				 propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
  
  for (int i = 0; i <= global.treemaxdepth; i ++) {
    s = s + mkformatstring("%7.3f ",
			   propose_depth[i] == 0 ? 0.0 : 100.0*(double)accept_depth[i]/(double)propose_depth[i]);
  }

  return s;
}

void
Death::initialize_mpi(MPI_Comm _communicator)
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
Death::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
Death::choose_death_location(int k,
			     double &ratio,
			     int &death_depth,
			     int &death_idx,
			     double &choose_prob,
			     int &death_valid)
{
  if (primary()) {
    //
    // Get tree structure ratio
    //
    
    if (hnk_get_kplus1_ratio(global.hnk, 
			     global.treemaxdepth, 
			     k - 1, 
			     &ratio) < 0) {
      ERROR("failed to get ratio for birth\n");
      return -1;
    }
    
    if (wavetree3d_sub_choose_death_global(global.wt, 
					   global.random.uniform(), 
					   global.treemaxdepth, 
					   &death_depth,
					   &death_idx,
					   &choose_prob) < 0) {
      death_valid = 0;
    }
  }
  
  return 0;
}

int
Death::communicate_death_location(int &death_valid,
				  int &death_idx,
				  int &death_depth)
{
  if (communicator != MPI_COMM_NULL) {

    if (MPI_Bcast(&death_valid, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast death valid\n");
    }
    
    if (death_valid) {
      if (MPI_Bcast(&death_idx, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast death index\n");
      }
      if (MPI_Bcast(&death_depth, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast death depth\n");
      }
    }
  }

  return 0;
}

int
Death::propose_death(int death_valid,
		     int death_idx,
		     int death_depth,
		     double &death_value)
{
  if (death_valid) {
    if (coefficient_histogram_propose_death(global.coeff_hist, death_idx) < 0) {
      ERROR("failed to update histogram for death proposal\n");
      return -1;
    }
    
    if (wavetree3d_sub_propose_death(global.wt,
				     death_idx, 
				     death_depth,
				     &death_value) < 0) {
      ERROR("failed to propose death\n");
      return -1;
    }
  }

  return 0;
}

int
Death::compute_reverse_death_probability(int death_idx,
					 int death_depth,
					 double death_value,
					 int &ii,
					 int &ij,
					 int &ik,
					 double &death_parent_coeff,
					 double &death_prob,
					 double &reverse_prob,
					 double &prior_prob)
{
  if (primary()) {
    
    if (wavetree3d_sub_3dindices(global.wt, death_idx, &ii, &ij, &ik) < 0) {
      ERROR("failed to compute 3d indices for death\n");
      return -1;
    }
    
    if (wavetree3d_sub_get_coeff(global.wt,
				 wavetree3d_sub_parent_index(global.wt, death_idx),
				 death_depth - 1,
				 &death_parent_coeff) < 0) {
      ERROR("failed to get parent coefficient for death\n");
      return -1;
    }
    
    if (wavetree_pp_death3d(global.proposal,
			    ii,
			    ij,
			    ik,
			    death_depth,
			    global.treemaxdepth,
			    death_parent_coeff,
			    death_value,
			    &death_prob) < 0) {
      ERROR("failed to get death probability\n");
      return -1;
    }
    
    if (wavetree3d_sub_reverse_death_global(global.wt, 
					    global.treemaxdepth,
					    death_depth,
					    death_idx, 
					    &reverse_prob) < 0) {
      ERROR("failed to reverse death (global)\n");
      return -1;
    }
    
    //
    // Compute the prior ratio
    //
    prior_prob = wavetree_pp_prior_probability3d(global.proposal,
						 ii, 
						 ij,
						 ik,
						 death_depth,
						 global.treemaxdepth,
						 death_parent_coeff,
						 death_value);
  }

  return 0;
}

int
Death::compute_likelihood(int death_idx,
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
Death::compute_acceptance(double proposed_likelihood,
			  double proposed_log_normalization,
			  double reverse_prob,
			  double choose_prob,
			  double death_prob,
			  double ratio,
			  double prior_prob,
			  bool &accept_proposal)
{
  if (primary()) {

    double u = log(global.random.uniform());
    
    accept_proposal = u < ((global.current_likelihood - proposed_likelihood)/global.temperature // Likelihood ratio
			   + log(reverse_prob) - 
			   log(choose_prob)                           // Node proposal ratio 
			   + log(death_prob)                          // Coefficient proposal ratio 
			   + log(ratio)                               // Tree prior ratio 
			   - log(prior_prob)                          // Coefficient prior ratio 
			   );

  }

  return 0;
}

int
Death::communicate_acceptance(bool &accept_proposal)
{
  if (communicator != MPI_COMM_NULL) {

    int ta;

    if (mpi_rank == 0) {
      ta = (int)accept_proposal;
    }

    if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast acceptance\n");
    }

    if (mpi_rank != 0) {
      accept_proposal = (bool)ta;
    }
  }

  return 0;
}
