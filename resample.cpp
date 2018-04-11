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

#include "resample.hpp"

extern "C" {
  #include "slog.h"
};

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

Resample::Resample(Global &_global) :
  global(_global),
  global_communicator(MPI_COMM_NULL),
  temperature_communicator(MPI_COMM_NULL),
  chain_communicator(MPI_COMM_NULL),
  resamplings(0),
  reselected(0),
  propagated(0),
  replaced(0),
  global_size(-1),
  global_rank(-1),
  temperature_size(-1),
  temperature_rank(-1),
  chain_size(-1),
  chain_rank(-1),
  processesperchain(-1),
  likelihood_array(nullptr),
  lambda_array(nullptr),
  weights_array(nullptr),
  sources_array(nullptr),
  requests_array(nullptr),
  send_buffer_size(-1),
  send_buffer(nullptr),
  recv_buffer_size(-1),
  recv_buffer(nullptr)
{
}

Resample::~Resample()
{
  delete [] send_buffer;
  delete [] recv_buffer;

  if (chain_rank == 0) {
    delete [] likelihood_array;
    delete [] weights_array;
    delete [] sources_array;
    delete [] requests_array;
  }
}

int
Resample::step(double resample_temperature)
{
  if (global_communicator == MPI_COMM_NULL) {
    ERROR("MPI Unitialized");
    return -1;
  }

  resamplings ++;
  int model_broadcast_required = 0;
  // double old_likelihood = global.current_likelihood;
  // double init_likelihood = global.likelihood_mpi();

  // INFO("Step %f %f", old_likelihood, init_likelihood);
  
  if (chain_rank == 0) {

    //
    // First distribute likelihoods between all chains
    //

    //
    // Strictly speaking, here we should include the normalizing constant of the hierarchical
    // parameters here as well but since we aren't doing hierarchical estimation here these
    // can be ommitted (all the same).
    //
    double like = global.current_likelihood;
    
    MPI_Allgather(&like, 1, MPI_DOUBLE, likelihood_array, 1, MPI_DOUBLE, temperature_communicator);

    double lambda = global.hierarchical->getparameter(0);

    MPI_Allgather(&lambda, 1, MPI_DOUBLE, lambda_array, 1, MPI_DOUBLE, temperature_communicator);

    // for (int i = 0; i < temperature_size; i ++) {
    //   INFO("%2d %f", i, likelihood_array[i]);
    // }
    
    //
    // Temper and normalize
    //
    double minnll = likelihood_array[0];
    for (int i = 1; i < temperature_size; i ++) {
      if (likelihood_array[i] < minnll) {
	minnll = likelihood_array[i];
      }
    }
    
    double sum = 0.0;
    for (int i = 0; i < temperature_size; i ++) {
      //
      // We subtract the minimum NLL and temper first before taking the exp
      // to try to prevent underflow.
      //
      weights_array[i] = (likelihood_array[i] - minnll)/resample_temperature;
      if (weights_array[i] > 700.0) {
	//
	// Above 700 we lose precision so we set probability to 0 here
	//
	weights_array[i] = 0.0;
      } else {
	weights_array[i] = exp(-weights_array[i]);
      }
      sum += weights_array[i];
    }

    for (int i = 0; i < temperature_size; i ++) {
      weights_array[i] /= sum;
    }

    //
    // Randomly select source of my resampled model
    //
    int new_source = global.random.select(temperature_size, weights_array);
    // INFO("Selected %d", new_source);

    //
    // Distribute sources
    //
    MPI_Allgather(&new_source, 1, MPI_INT, sources_array, 1, MPI_INT, temperature_communicator);

    //
    // Encode my model into the send buffer
    //
    int send_length = wavetree3d_sub_encode(global.wt, send_buffer, send_buffer_size);
    if (send_length < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to encode wavetree\n");
    }
    
    //
    // Loop through source array and send my model to those who need it
    //
    int nsends = 0;
    for (int i = 0; i < temperature_size; i ++) {
      if (i != temperature_rank && sources_array[i] == temperature_rank) {

	// INFO("Sending to %d", i);
	MPI_Isend(send_buffer, send_length, MPI_BYTE, i, 0, temperature_communicator, &requests_array[nsends]);
	nsends ++;
	propagated ++;

      }
    }
    
    //
    // Receive my model (if I didn't select myself)
    //
    if (sources_array[temperature_rank] != temperature_rank) {

      replaced ++;
      
      MPI_Status status;
      MPI_Probe(sources_array[temperature_rank], 0, temperature_communicator, &status);

      int recv_length;
      MPI_Get_count(&status, MPI_BYTE, &recv_length);

      if (recv_length > recv_buffer_size) {
	throw WAVETOMO2DEXCEPTION("Too many bytes for buffer: %d > %d", recv_length, recv_buffer_size);
      }
      
      MPI_Recv(recv_buffer, recv_length, MPI_BYTE, sources_array[temperature_rank], 0, temperature_communicator, &status);

      //
      // Set my model and the current likelihood
      //
      if (wavetree3d_sub_decode(global.wt, recv_buffer, recv_length) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to decode wavetree\n");
      }
      global.current_likelihood = likelihood_array[sources_array[temperature_rank]];
      global.hierarchical->setparameter(0, lambda_array[sources_array[temperature_rank]]);

      model_broadcast_required = 1;
    } else {
      //
      // Kept current model
      //
      reselected ++;
    }

    //
    // Wait for all non-blocking sends to finish
    //
    for (int i = 0; i < nsends; i ++) {
      MPI_Status status;
      
      MPI_Wait(&requests_array[i], &status);
    }

    if (model_broadcast_required) {
      // INFO("New Likelihood: %f", global.current_likelihood);
    }

  }

  if (processesperchain > 1) {

    //
    // Broadcast new model to rest of chain parallel processes if required
    //
    MPI_Bcast(&model_broadcast_required, 1, MPI_INT, 0, chain_communicator);

    if (model_broadcast_required) {
      double lambda_scale;
      
      MPI_Bcast(&global.current_likelihood, 1, MPI_DOUBLE, 0, chain_communicator);
      
      if (chain_rank == 0) {

	lambda_scale = global.hierarchical->getparameter(0);
	MPI_Bcast(&lambda_scale, 1, MPI_DOUBLE, 0, chain_communicator);

	//
	// Reencode my model and broadcast
	//
	int send_length = wavetree3d_sub_encode(global.wt, send_buffer, send_buffer_size);
	if (send_length < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to encode wavetree\n");
	}

	MPI_Bcast(&send_length, 1, MPI_INT, 0, chain_communicator);
	MPI_Bcast(send_buffer, send_length, MPI_BYTE, 0, chain_communicator);

	
      } else {

	MPI_Bcast(&lambda_scale, 1, MPI_DOUBLE, 0, chain_communicator);
	global.hierarchical->setparameter(0, lambda_scale);
	
	int recv_length;
	
	MPI_Bcast(&recv_length, 1, MPI_INT, 0, chain_communicator);
	MPI_Bcast(recv_buffer, recv_length, MPI_BYTE, 0, chain_communicator);

	if (wavetree3d_sub_decode(global.wt, recv_buffer, recv_length) < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to decode wavetree\n");
	}
      }

    }

  }

  // double new_likelihood = global.likelihood_mpi();

  // INFO("Old Like %f Current %f New Like %f", old_likelihood, global.current_likelihood, new_likelihood);
  
  // if (fabs(new_likelihood - global.current_likelihood) > 1.0e-9) {
  //   throw WAVETOMO2DEXCEPTION("Likelihood mismatch %.9g", fabs(global.current_likelihood - new_likelihood));
  // }

  return model_broadcast_required;
}

std::string
Resample::write_short_stats()
{
  return mkformatstring("ReSmpl %6d: %6d K %6d P %6d R",
			resamplings,
			reselected,
			propagated,
			replaced);
}

std::string
Resample::write_long_stats()
{
  return write_short_stats();
}

void
Resample::initialize_mpi(MPI_Comm _global_communicator,
			 MPI_Comm _temperature_communicator,
			 MPI_Comm _chain_communicator)
{
  if (MPI_Comm_dup(_global_communicator, &global_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to duplicate global communicator");
  }
  
  if (MPI_Comm_size(global_communicator, &global_size) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure");
  }
  if (MPI_Comm_rank(global_communicator, &global_rank) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure");
  }
  
  if (MPI_Comm_dup(_temperature_communicator, &temperature_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to duplicate global communicator");
  }
  if (MPI_Comm_dup(_chain_communicator, &chain_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to duplicate global communicator");
  }

  if (MPI_Comm_rank(chain_communicator, &chain_rank) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure");
  }
    
  if (chain_rank == 0) {
    if (MPI_Comm_rank(temperature_communicator, &temperature_rank) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("MPI Failure");
    }
      
    if (MPI_Comm_size(temperature_communicator, &temperature_size) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("MPI Failure");
    }
      
    likelihood_array = new double[temperature_size];
    lambda_array = new double[temperature_size];
    weights_array = new double[temperature_size];
    sources_array = new int[temperature_size];
    requests_array = new MPI_Request[temperature_size];

    processesperchain = global_size/temperature_size;
    
  }

  if (MPI_Bcast(&temperature_size, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to broadcast ntotalchains.");
  }
  if (MPI_Bcast(&processesperchain, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to broardcast processesperchain.");
  }

  INFO("%03d: Global Size: %d NChains: %d PPC: %d\n", global_rank, global_size, temperature_size, processesperchain);

  send_buffer_size = global.ncoeff * 3 * sizeof(double);
  recv_buffer_size = send_buffer_size;

  send_buffer = new char[send_buffer_size];
  recv_buffer = new char[recv_buffer_size];
}


