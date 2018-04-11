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

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "ptexchangeslice.hpp"

PTExchangeSlice::PTExchangeSlice(GlobalSlice &_global) :
  global(_global),
  propose(0),
  accept(0),
  global_communicator(MPI_COMM_NULL),
  global_size(-1),
  global_rank(-1),
  temperature_communicator(MPI_COMM_NULL),
  temperature_size(-1),
  temperature_rank(-1),
  chain_communicator(MPI_COMM_NULL),
  chain_rank(-1),
  ntotalchains(-1),
  processesperchain(-1),
  ptpairs(nullptr),
  partner(-1),
  send(false),
  ptaccept(false),
  send_length(-1),
  recv_length(-1),
  partner_likelihood(-1.0),
  partner_temperature(-1.0),
  u(-1.0),
  exchanged(-1),
  send_buffer_size(-1),
  send_buffer(nullptr),
  recv_buffer_size(-1),
  recv_buffer(nullptr)
{
  for (int i = 0; i < 3; i ++) {
    sendmsg[i] = 0.0;
    recvmsg[i] = 0.0;
  }
}

PTExchangeSlice::~PTExchangeSlice()
{
  delete [] ptpairs;
  
  delete [] send_buffer;
  delete [] recv_buffer;
}

int
PTExchangeSlice::step()
{
  if (global_communicator == MPI_COMM_NULL) {
    return -1;
  }

  // double old_likelihood = global.current_likelihood;
  // double init_likelihood = global.likelihood_mpi();
  // INFO("Step %f %f", old_likelihood, init_likelihood);
  
  //
  // Generate the list of process pairs to generate exchanges
  //
  ptaccept = false;
  if (chain_rank == 0) {
    propose ++;
    
    if (temperature_rank == 0) {

      //
      // We do (potentially) 2 lots of shuffling here to ensure that we attempt an exchange
      // between chains at different temperatures rather than wasting exchanges between
      // the same temperature.
      //
      for (int j = 0; j < temperature_size; j ++) {
	//
	// Place all the mpi_rank's of chain_rank == 0 processes in the list.
	//
	transposed_ptpairs[j] = j;
      }

      //
      // Shuffle in temperature if necessary
      //
      if (chainspertemperature > 1) {
	for (int j = 0; j < ntemperatures; j ++) {
	  global.random.shuffle(chainspertemperature, transposed_ptpairs + j*chainspertemperature);
	}
      }

      //
      // Transpose
      //
      for (int i = 0; i < ntemperatures; i ++) {
	for (int j = 0; j < chainspertemperature; j ++) {

	  ptpairs[j * ntemperatures + i] = transposed_ptpairs[i * chainspertemperature + j];

	}
      }

      //
      // Shuffle between temperatures
      //
      for (int j = 0; j < chainspertemperature; j ++) {
	global.random.shuffle(ntemperatures, ptpairs + j * ntemperatures);
      }
      
    }

    if (MPI_Bcast(ptpairs, temperature_size, MPI_INT, 0, temperature_communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broadcast ptpairs\n");
    }

    //
    // Find self in shuffled list
    //
    partner = -1;
    for (int j = 0; j < temperature_size; j ++) {
      if (ptpairs[j] == temperature_rank) {
	if (j % 2 == 0) {
	  send = true;
	  partner = ptpairs[j + 1];
	} else {
	  send = false;
	  partner = ptpairs[j - 1];
	}
	break;
      }
    }
    
    if (partner < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to find self in exchange list\n");
    }
    
    //
    // Exchange information with exchange partner
    //
    current_logprior = wavetree2d_sub_logpriorprobability(global.wt, global.proposal);
    
    sendmsg[0] = global.current_likelihood;
    sendmsg[1] = global.current_log_normalization;
    sendmsg[2] = global.temperature;
    sendmsg[3] = current_logprior;

    if (send) {
      
      sendmsg[4] = global.random.uniform();
      
      if (MPI_Send(sendmsg, 5, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
      	throw WAVETOMO2DEXCEPTION("Failed to send current state to partner seconday\n");
      }      

      if (MPI_Recv(recvmsg, 5, MPI_DOUBLE, partner, 0, temperature_communicator, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      	throw WAVETOMO2DEXCEPTION("Failed to recv current state to partner secondary\n");
      }	
      
      partner_likelihood = recvmsg[0];
      partner_log_normalization = recvmsg[1];
      partner_temperature = recvmsg[2];
      partner_logprior = recvmsg[3];
      u = sendmsg[4];
    } else {

      sendmsg[4] = 0.0;

      if (MPI_Recv(recvmsg, 5, MPI_DOUBLE, partner, 0, temperature_communicator, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      	throw WAVETOMO2DEXCEPTION("Failed to recv current state to partner primary\n");
      }

      if (MPI_Send(sendmsg, 5, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
      	throw WAVETOMO2DEXCEPTION("Failed to send current state to partner primary\n");
      }      
      
      partner_likelihood = recvmsg[0];
      partner_log_normalization = recvmsg[1];
      partner_temperature = recvmsg[2];
      partner_logprior = recvmsg[3];
      u = recvmsg[4];
      
    }
    
    ptaccept = log(u) < ((((current_logprior + global.current_likelihood + global.current_log_normalization) -
			   (partner_logprior + partner_likelihood + partner_log_normalization))/global.temperature) +
			 
     			 (((partner_logprior + partner_likelihood + partner_log_normalization) -
			   (current_logprior + global.current_likelihood + global.current_log_normalization))/partner_temperature));
    
    if (ptaccept) {

      accept ++;
      
      //
      // Send model
      //
      if (send) {
	
	send_length = wavetree2d_sub_encode(global.wt, send_buffer, send_buffer_size);
	
	if (send_length < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to encode wavetree\n");
	}
	
	if (MPI_Send(&send_length, 1, MPI_INT, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send length to partner secondary\n");
	}
	if (MPI_Send(send_buffer, send_length, MPI_BYTE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send encoded model to partner seconday\n");
	}
	
	if (MPI_Recv(&recv_length, 1, MPI_INT, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv length to partner secondary\n");
	}
	if (MPI_Recv(recv_buffer, recv_length, MPI_BYTE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv encoded model to partner seconday\n");
	}
	
      } else {
	
	if (MPI_Recv(&recv_length, 1, MPI_INT, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv length to partner primary\n");
	}
	if (MPI_Recv(recv_buffer, recv_length, MPI_BYTE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv encoded model to partner primary\n");
	}
	  
	
	send_length = wavetree2d_sub_encode(global.wt, send_buffer, send_buffer_size);
	
	if (send_length < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to encode wavetree\n");
	}
	
	if (MPI_Send(&send_length, 1, MPI_INT, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send length to partner secondary\n");
	}
	if (MPI_Send(send_buffer, send_length, MPI_BYTE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send encoded model to partner seconday\n");
	}
	
      }
      
      if (wavetree2d_sub_decode(global.wt, recv_buffer, recv_length) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to decode wavetree\n");
      }

      //
      // Send hierarchical scale
      //
      double new_lambda_scale = 0.0;
      
      if (send) {

	double lambda_scale = global.hierarchical->getparameter(0);
	
	if (MPI_Send(&lambda_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send lambda\n");
	}
	
	if (MPI_Recv(&new_lambda_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv lambda\n");
	}
	
      } else {
	
	double lambda_scale = global.hierarchical->getparameter(0);

	if (MPI_Recv(&new_lambda_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv lambda\n");
	}
	  
	if (MPI_Send(&lambda_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send length to partner secondary\n");
	}
      }

      //
      // Send hierarchical prior scale
      //
      double new_prior_scale = 0.0;
      
      if (send) {

	double prior_scale;
	wavetree_pp_setscale(global.proposal, 0.0, &prior_scale);
	
	if (MPI_Send(&prior_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send prior\n");
	}
	
	if (MPI_Recv(&new_prior_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv prior\n");
	}
	
      } else {
	
	double prior_scale;
	wavetree_pp_setscale(global.proposal, 0.0, &prior_scale);

	if (MPI_Recv(&new_prior_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator, &status) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to recv prior\n");
	}
	  
	if (MPI_Send(&prior_scale, 1, MPI_DOUBLE, partner, 0, temperature_communicator) != MPI_SUCCESS) {
	  throw WAVETOMO2DEXCEPTION("Failed to send length to partner secondary\n");
	}
      }

      global.hierarchical->setparameter(0, new_lambda_scale);
      global.current_likelihood = partner_likelihood;
      global.current_log_normalization = partner_log_normalization;

      if (wavetree_pp_setscale(global.proposal, new_prior_scale, NULL) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to set new prior scale\n");
      }
    }
  }
  
  //
  // Now we may need to distribute exhanged models to each process in our chain
  //
  exchanged = 0;
  if (chain_rank == 0 && ptaccept) {
    exchanged = 1;
  }
  
  if (processesperchain > 1) {
    
    if (MPI_Bcast(&exchanged, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to broad cast exchanged\n");
    }
    
    if (exchanged) {

      double lambda_scale;
      double prior_scale;
      
      if (chain_rank == 0) {
	lambda_scale = global.hierarchical->getparameter(0);
	wavetree_pp_setscale(global.proposal, 0.0, &prior_scale);
	
	send_length = wavetree2d_sub_encode(global.wt, send_buffer, send_buffer_size);
	if (send_length < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to encode wavetree\n");
	}
      }

      if (MPI_Bcast(&global.current_likelihood, 1, MPI_DOUBLE, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad cast likelihood\n");
      }

      if (MPI_Bcast(&global.current_log_normalization, 1, MPI_DOUBLE, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad cast log normalization\n");
      }

      if (MPI_Bcast(&lambda_scale, 1, MPI_DOUBLE, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad cast lambda\n");
      }
      
      if (MPI_Bcast(&prior_scale, 1, MPI_DOUBLE, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad cast lambda\n");
      }

      if (MPI_Bcast(&send_length, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad cast send length\n");
      }
      
      if (MPI_Bcast(send_buffer, send_length, MPI_BYTE, 0, chain_communicator) != MPI_SUCCESS) {
	throw WAVETOMO2DEXCEPTION("Failed to broad case encoded model\n");
      }
      
      if (chain_rank != 0) {

	global.hierarchical->setparameter(0, lambda_scale);

	if (wavetree_pp_setscale(global.proposal, prior_scale, NULL) < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to set prior scale\n");
	}
	
	if (wavetree2d_sub_decode(global.wt, send_buffer, send_length) < 0) {
	  throw WAVETOMO2DEXCEPTION("Failed to decode wavetree\n");
	}
      }
    }
  }

  // double new_likelihood = global.likelihood_mpi();

  // INFO("%d Old Like %f Current %f New Like %f", exchanged, old_likelihood, global.current_likelihood, new_likelihood);
  
  // if (new_likelihood != global.current_likelihood) {
  //   throw WAVETOMO2DEXCEPTION("Likelihood mismatch %.9g", fabs(global.current_likelihood - new_likelihood));
  // }

  return exchanged;
}


std::string
PTExchangeSlice::write_short_stats()
{
  return mkformatstring("PTxchg %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}


std::string
PTExchangeSlice::write_long_stats()
{
  return write_short_stats();
}

void
PTExchangeSlice::initialize_mpi(MPI_Comm _global_communicator,
				MPI_Comm _temperature_communicator,
				MPI_Comm _chain_communicator,
				int _ntemperatures)
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
      
    ntotalchains = temperature_size;
    if (global_size % ntotalchains != 0) {
      throw WAVETOMO2DEXCEPTION("Global size/Total Chains mismatch: %d %d\n", global_size, ntotalchains);
    }
    
    processesperchain = global_size/ntotalchains;
    chainspertemperature = ntotalchains/_ntemperatures;
    ntemperatures = _ntemperatures;
    
    ptpairs = new int[temperature_size];
    transposed_ptpairs = new int[temperature_size];
  }

  
  if (MPI_Bcast(&ntotalchains, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to broadcast ntotalchains.");
  }
  if (MPI_Bcast(&processesperchain, 1, MPI_INT, 0, chain_communicator) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("Failed to broardcast processesperchain.");
  }

  INFO("%03d: Global Size: %d NChains: %d PPC: %d\n", global_rank, global_size, ntotalchains, processesperchain);

  send_buffer_size = global.ncoeff * 3 * sizeof(double);
  recv_buffer_size = send_buffer_size;

  send_buffer = new char[send_buffer_size];
  recv_buffer = new char[recv_buffer_size];

  
}

  
