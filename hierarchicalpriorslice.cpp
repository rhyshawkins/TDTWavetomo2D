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

extern "C" {
#include "slog.h"
};

#include "hierarchicalpriorslice.hpp"

#include "wavetomo2dexception.hpp"
#include "wavetomo2dutil.hpp"

HierarchicalPriorSlice::HierarchicalPriorSlice(GlobalSlice &_global, double _sigma) :
  global(_global),
  sigma(_sigma),
  propose(0),
  accept(0),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1),
  alpha(1.0),
  beta(0.0)
{
  double initialb;
  wavetree_pp_setscale(_global.proposal, 0.0, &initialb);

  //
  // Set the mode of the inverse gamma distribution to be the initial value
  //
  beta = (alpha + 1.0) * initialb;
}

HierarchicalPriorSlice::~HierarchicalPriorSlice()
{
}

int
HierarchicalPriorSlice::step()
{
  propose ++;

  double value;
  double value_prior_ratio = 0.0;
  int valid_proposal = 0;
  double old_value;

  //
  // Retrieve the old value
  //
  wavetree_pp_setscale(global.proposal, 0.0, &old_value);
  value = old_value;

  if (choose_value(value,
		   value_prior_ratio,
		   valid_proposal) < 0) {
    return -1;
  }
  
  if (communicate_value(valid_proposal,
			value) < 0) {
    return -1;
  }

  memset(&last_step, 0, sizeof(chain_history_change_t));

  //
  // hierarchical prior scaling proposal are placed into hyper steps
  //
  last_step.header.type = CH_HYPER;
  
  last_step.perturbation.hyper.index = 0;
  last_step.perturbation.hyper.new_value = value;
  last_step.perturbation.hyper.old_value = old_value;
  
  if (valid_proposal) {

    double current_log_prior = wavetree2d_sub_logpriorprobability(global.wt, global.proposal);
    if (wavetree_pp_setscale(global.proposal, value, NULL) < 0) {
      ERROR("Failed to set prior scale");
      return -1;
    }

    double proposed_log_prior = wavetree2d_sub_logpriorprobability(global.wt, global.proposal);
    if (wavetree_pp_setscale(global.proposal, old_value, NULL) < 0) {
      ERROR("Failed to reset prior scale");
      return -1;
    }

    
    bool accept_proposal = false;
    if (compute_acceptance(value_prior_ratio,
			   current_log_prior,
			   proposed_log_prior,
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

      last_step.header.accepted = 1;
      
      if (wavetree_pp_setscale(global.proposal, value, NULL) < 0) {
	ERROR("Failed to set prior scale");
	return -1;
      }
      
      return 1;
      
    } else {
      
      //
      // Reject
      //

      return 0;
    }

  }
  
  return 0;
}

std::string
HierarchicalPriorSlice::write_short_stats()
{
  return mkformatstring("HierarchicalPriorSlice %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
HierarchicalPriorSlice::write_long_stats()
{
  return write_short_stats();
}

void
HierarchicalPriorSlice::initialize_mpi(MPI_Comm _communicator)
{
  MPI_Comm_dup(_communicator, &communicator);

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }
}

void
HierarchicalPriorSlice::get_last_step(chain_history_change_t *_last_step)
{
  memcpy(_last_step, &last_step, sizeof(chain_history_change_t));
}

bool
HierarchicalPriorSlice::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
HierarchicalPriorSlice::choose_value(double &value,
				     double &log_value_prior_ratio,
				     int &valid_proposal)
{
  if (primary()) {
    double old_value = value;
    value *= exp(sqrt(global.temperature) * global.random.normal(sigma));
    
    valid_proposal = 1;
    log_value_prior_ratio = beta * (old_value*old_value - value*value);
  }

  return 0;
}
  

int
HierarchicalPriorSlice::communicate_value(int &valid_proposal,
					  double &value)
{
  if (communicator != MPI_COMM_NULL) {

    if (MPI_Bcast(&valid_proposal, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast valid proposal\n");
    }
    
    if (valid_proposal) {

      if (MPI_Bcast(&value, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broadcast value\n");
      }
      
    }
  }

  return 0;
}
  
int
HierarchicalPriorSlice::compute_acceptance(double log_value_prior_ratio,
					   double current_log_prior,
					   double proposed_log_prior,
					   bool &accept_proposal)
{
  if (primary()) {
    
    double u = log(global.random.uniform());

    //
    // There is an option to temper the prior using global.temperature here but
    // Jan suggests exposing the full prior is better for PT
    //
    double alpha =
      log_value_prior_ratio +
      proposed_log_prior - current_log_prior;

    accept_proposal = (u < alpha);

  }

  return 0;
}

int
HierarchicalPriorSlice::communicate_acceptance(bool &accept_proposal)
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
