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

#include <getopt.h>

#include <gmp.h>

#include <mpi.h>

extern "C" {
#include "wavetree3d_sub.h"
#include "wavetreepp.h"
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "slog.h"
};

#include "global.hpp"
#include "birth.hpp"
#include "death.hpp"
#include "value.hpp"
#include "hierarchical.hpp"
#include "hierarchicalprior.hpp"
#include "ptexchange.hpp"

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

static char short_options[] = "i:I:M:H:L:o:x:y:z:u:n:N:a:A:t:S:F:l:k:B:Pw:W:v:c:T:m:e:rU:R:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"prior-file", required_argument, 0, 'M'},
  {"hierarchical-sigma", required_argument, 0, 'H'},
  {"hierarchical-prior-sigma", required_argument, 0, 'L'},
  {"output", required_argument, 0, 'o'},
  
  {"degree-x", required_argument, 0, 'x'},
  {"degree-y", required_argument, 0, 'y'},
  {"degree-z", required_argument, 0, 'z'},
  {"super-resolution", required_argument, 0, 'u'},

  {"lonmin", required_argument, 0, 'n'},
  {"lonmax", required_argument, 0, 'N'},
  {"latmin", required_argument, 0, 'a'},
  {"latmax", required_argument, 0, 'A'},

  {"total", required_argument, 0, 't'},
  {"seed", required_argument, 0, 'S'},
  {"seed-mult", required_argument, 0, 'F'},

  {"lambda", required_argument, 0, 'l'},

  {"kmax", required_argument, 0, 'k'},

  {"birth-probability", required_argument, 0, 'B'},

  {"posteriork", no_argument, 0, 'P'},

  {"wavelet-xy", required_argument, 0, 'w'},
  {"wavelet-z", required_argument, 0, 'W'},

  {"verbosity", required_argument, 0, 'v'},

  {"chains", required_argument, 0, 'c'},
  {"temperatures", required_argument, 0, 'T'},
  {"max-temperature", required_argument, 0, 'm'},

  {"exchange-rate", required_argument, 0, 'e'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static constexpr double PB_MAX = 0.45;

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Parameters
  //
  char *input_obs;
  char *initial_model;
  char *prior_file;
  double hierarchical_sigma;
  double hierarchical_prior_sigma;
  char *output_prefix;

  int degreex;
  int degreey;
  int degreez;
  int super_resolution;

  double lonmin;
  double lonmax;
  double latmin;
  double latmax;

  int total;
  int seed_base;
  int seed_mult;

  double lambda;
  int kmax;

  double Pb;

  bool posteriork;

  int wavelet_xy;
  int wavelet_z;

  int verbosity;

  int chains;
  int temperatures;
  int exchange_rate;
  double max_temperature;

  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  //
  // Defaults
  //

  input_obs = nullptr;
  initial_model = nullptr;
  prior_file = nullptr;
  hierarchical_sigma = 0.0;
  hierarchical_prior_sigma = 0.0;
  output_prefix = nullptr;

  degreex = 7;
  degreey = 6;
  degreez = 5;
  super_resolution = 0;
  
  lonmin = -10.0;
  lonmax = 10.0;
  latmin = -10.0;
  latmax = 10.0;

  total = 10000;
  seed_base = 983;
  seed_mult = 101;

  lambda = 1.0;
  kmax = 100;

  Pb = 0.05;

  posteriork = false;

  wavelet_xy = 0;
  wavelet_z = 0;

  verbosity = 1000;

  chains = 1;
  temperatures = 1;
  max_temperature = 1000.0;
  exchange_rate = 10;

  //
  // Command line parameters
  //
  option_index = 0;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input_obs = optarg;
      break;

    case 'I':
      initial_model = optarg;
      break;

    case 'M':
      prior_file = optarg;
      break;

    case 'H':
      hierarchical_sigma = atof(optarg);
      if (hierarchical_sigma < 0.0) {
	fprintf(stderr, "error: hierarchical sigma must be 0 (hierarchical sampling disabled) or greater\n");
	return -1;
      }
      break;

    case 'L':
      hierarchical_prior_sigma = atof(optarg);
      if (hierarchical_prior_sigma < 0.0) {
	fprintf(stderr, "error: hierarchical prior sigma must be greater than 0\n");
	return -1;
      }
      break;

    case 'o':
      output_prefix = optarg;
      break;

    case 'x':
      degreex = atoi(optarg);
      if (degreex < 1 || degreex > 16) {
	fprintf(stderr, "error: degree x must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'y':
      degreey = atoi(optarg);
      if (degreey < 1 || degreey > 16) {
	fprintf(stderr, "error: degree y must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'z':
      degreez = atoi(optarg);
      if (degreez < 1 || degreez > 16) {
	fprintf(stderr, "error: degree z must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'u':
      super_resolution = atoi(optarg);
      if (super_resolution < 0) {
	fprintf(stderr, "error: super resolution must be 0 or greater\n");
	return -1;
      }
      break;
      
    case 'n':
      lonmin = atof(optarg);
      break;

    case 'N':
      lonmax = atof(optarg);
      break;

    case 'a':
      latmin = atof(optarg);
      break;

    case 'A':
      latmax = atof(optarg);
      break;

    case 't':
      total = atoi(optarg);
      if (total < 1) {
	fprintf(stderr, "error: total must be greater than 0\n");
	return -1;
      }
      break;

    case 'S':
      seed_base = atoi(optarg);
      break;

    case 'F':
      seed_mult = atoi(optarg);
      if (seed_mult <= 0) {
	fprintf(stderr, "error: seed multiplier must be greater than 0\n");
	return -1;
      }
      break;

    case 'l':
      lambda = atof(optarg);
      if (lambda <= 0.0) {
	fprintf(stderr, "Lambda std dev must be greater than 0\n");
	return -1;
      }
      break;
      
    case 'k':
      kmax = atoi(optarg);
      if (kmax < 1) {
	fprintf(stderr, "error: kmax must be greater than 0\n");
	return -1;
      }
      break;

    case 'B':
      Pb = atof(optarg);
      if (Pb < 0.0 || Pb > PB_MAX) {
	fprintf(stderr, "error: birth probability must be between 0 and %.3f\n", PB_MAX);
	return -1;
      }
      break;

    case 'P':
      posteriork = true;
      break;

    case 'w':
      wavelet_xy = atoi(optarg);
      if (wavelet_xy < 0 || wavelet_xy > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'W':
      wavelet_z = atoi(optarg);
      if (wavelet_z < 0 || wavelet_z > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'v':
      verbosity = atoi(optarg);
      break;

    case 'c':
      chains = atoi(optarg);
      if (chains < 1) {
	fprintf(stderr, "error: no. chains must be greater than 0\n");
	return -1;
      }
      break;

    case 'T':
      temperatures = atoi(optarg);
      if (temperatures < 1) {
	fprintf(stderr, "error: no. temperatures must be greater than 0\n");
	return -1;
      }
      break;

    case 'm':
      max_temperature = atof(optarg);
      if (max_temperature < 1.0) {
	fprintf(stderr, "error: maximum temperature must be 1.0 or greater\n");
	return -1;
      }
      break;

    case 'e':
      exchange_rate = atoi(optarg);
      if (exchange_rate < 0) {
	fprintf(stderr, "error: exchange rate must be 0 or greater\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  std::string logfile = mkfilenamerank(output_prefix, "log.txt", mpi_rank);
  if (slog_set_output_file(logfile.c_str(),
			   SLOG_FLAGS_CLEAR) < 0) {
    fprintf(stderr, "error: failed to redirect log file\n");
    return -1;
  }

  if (input_obs == nullptr) {
    ERROR("error: required input parameter input observations missing\n");
    return -1;
  }

  if (prior_file == nullptr) {
    ERROR("error: required prior file parameter missing\n");
    return -1;
  }

  //
  // The chains variable specifies the number of chains per temperature so that we
  // have temperatures * chains total chains. This must be a factor of the mpi size
  // with the factorization being temperatures * chains * processesperchain == mpi_size.
  //
  // For example, with mpi_size = 16, chains = 4, temperatures = 2, processesperchain
  // would be 2 and the layout would be:
  //
  // mpi_rank   chainid  chainrank  temperatureid
  //        0         0          0              0
  //        1         0          1              0
  //        2         1          0              0
  //        3         1          1              0
  //        4         2          0              0
  //        5         2          1              0
  //        6         3          0              0
  //        7         3          1              0
  //        8         4          0              1
  //        9         4          1              1
  //       10         5          0              1
  //       11         5          1              1
  //       12         6          0              1
  //       13         6          1              1
  //       14         7          0              1
  //       15         7          1              1
  //


  int ntotalchains = temperatures * chains;
  if (ntotalchains == 0 || mpi_size % ntotalchains != 0) {
    ERROR("error: no. temperatures and no. chains incompatible with mpi size: %d x %d = %d : %d\n",
	    temperatures, chains, ntotalchains, mpi_size);
    return -1;
  }

  if (temperatures > 1 && ntotalchains % 2 != 0) {
    ERROR("error: no. total chains (no. temperatures * no. chains) must be even\n");
    return -1;
  }
  
  int processesperchain = mpi_size/ntotalchains;
  int chain_id = mpi_rank/processesperchain;
  int chain_rank = mpi_rank % processesperchain;
  int temp_id = chain_id/chains;
  double temperature;
  if (temperatures == 1) {
    temperature = 1.0;
  } else {
    temperature = pow(10.0, log10(max_temperature) * (double)temp_id/(double)(temperatures - 1));
  }

  const char *initial_model_ptr = nullptr;
  std::string initial_model_rank;
  if (initial_model != nullptr) {
    initial_model_rank = mkfilenamerank(initial_model, "final_model.txt", chain_id);
    initial_model_ptr = initial_model_rank.c_str();
  }

  Global *global = new Global(input_obs,
			      initial_model_ptr,
			      prior_file,
			      degreex,
			      degreey,
			      degreez,
			      super_resolution,
			      lonmin,
			      lonmax,
			      latmin,
			      latmax,
			      seed_base + mpi_rank*seed_mult,
			      kmax,
			      lambda,
			      posteriork,
			      wavelet_xy,
			      wavelet_z);

  Birth *birth = new Birth(*global);
  Death *death = new Death(*global);
  Value *value = new Value(*global);

  Hierarchical *hierarchical = nullptr;
  if (hierarchical_sigma > 0.0) {
    hierarchical = new Hierarchical(*global, hierarchical_sigma);
  }

  HierarchicalPrior *hierarchical_prior = nullptr;
  if (hierarchical_prior_sigma > 0.0) {
    hierarchical_prior = new HierarchicalPrior(*global, hierarchical_prior_sigma);
  }

  PTExchange *ptexchange = nullptr;
  if (temperatures > 1) {
    ptexchange = new PTExchange(*global);
  }
  
  MPI_Comm chain_communicator;
  MPI_Comm temperature_communicator;

  //
  // Chain communicator used for parallel likelihood evaluation
  //
  MPI_Comm_split(MPI_COMM_WORLD, chain_id, mpi_rank, &chain_communicator);
  MPI_Comm_set_errhandler(chain_communicator, MPI_ERRORS_RETURN);

  global->initialize_mpi(chain_communicator, temperature);
  birth->initialize_mpi(chain_communicator);
  death->initialize_mpi(chain_communicator);
  value->initialize_mpi(chain_communicator);
  if (hierarchical) {
    hierarchical->initialize_mpi(chain_communicator);
  }

  //
  // Temperature communicator used for orchestrating PT exchanges
  //
  if (ptexchange != nullptr) {
    MPI_Comm_split(MPI_COMM_WORLD, chain_rank == 0, mpi_rank, &temperature_communicator);
    MPI_Comm_set_errhandler(temperature_communicator, MPI_ERRORS_RETURN);
    int temperature_rank;
    MPI_Comm_rank(temperature_communicator, &temperature_rank);
    if (mpi_rank == 0) {
      if (temperature_rank != 0) {
	throw WAVETOMO2DEXCEPTION("MPI Rank unexpected: %d != %d\n", mpi_rank, temperature_rank);
      }
    }
    
    ptexchange->initialize_mpi(MPI_COMM_WORLD,
			       temperature_communicator,
			       chain_communicator,
			       temperatures);
  }
  
  global->current_likelihood = global->likelihood_mpi(global->current_log_normalization);
    
  if (chain_rank == 0) {
    INFO("%03d Initial Likelihood: %f (%f)\n", chain_id, global->current_likelihood, global->current_log_normalization);
  }

  //
  // This initializes the residual tracking code.
  //
  global->accept();

  int *khistogram = nullptr;
  if (chain_rank == 0) {
    khistogram = new int[kmax];
    for (int i = 0; i < kmax; i ++) {
      khistogram[i] = 0;
    }
  }

  FILE *fp_ch = NULL;
  if (!posteriork && chain_rank == 0) {
    if (chain_history_initialise(global->ch,
				 wavetree3d_sub_get_S_v(global->wt),
				 global->current_likelihood,
				 global->temperature,
				 global->hierarchical->getparameter(0)) < 0) {
      ERROR("error: failed to initialise chain history\n");
      return -1;
    }

    std::string filename = mkfilenamerank(output_prefix, "ch.dat", chain_id);
    fp_ch = fopen(filename.c_str(), "w");
    if (fp_ch == NULL) {
      ERROR("error: failed to create chain history file\n");
      return -1;
    }
  }

  if (chain_rank == 0) {
    INFO("Starting Iterations");
  }
  
  for (int i = 0; i < total; i ++) {

    //
    // Make sure we're all here
    //
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Failed to barrier\n");
    }
  
    double u;
    if (chain_rank == 0) {
      u = global->random.uniform();
    }

    MPI_Bcast(&u, 1, MPI_DOUBLE, 0, chain_communicator);
    
    if (u < Pb) {

      //
      // Birth
      //
      if (birth->step() < 0) {
	ERROR("error: failed to do birth step\n");
	return -1;
      }

    } else if (u < (2.0 * Pb)) {

      //
      // Death
      //
      if (death->step() < 0) {
	ERROR("error: failed to do death step\n");
	return -1;
      }

    } else {

      //
      // Value
      //
      
      if (value->step() < 0) {
	ERROR("error: failed to do value step\n");
	return -1;
      }

    }

    int current_k = wavetree3d_sub_coeff_count(global->wt);

    if (chain_rank == 0) {
      khistogram[current_k - 1] ++;

      if (!posteriork) {
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}

	chain_history_change_t step;
      
	if (wavetree3d_sub_get_last_perturbation(global->wt, &step) < 0) {
	  ERROR("error: failed to get last step\n");
	  return -1;
	}
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = 1.0;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }

    //
    // Hierarchical
    //
    if (hierarchical != nullptr) {
      if (hierarchical->step() < 0) {
	ERROR("error: failed to do hierarchical step");
	return -1;
      }

      if (chain_rank == 0) {
	chain_history_change_t step;
	
	hierarchical->get_last_step(&step);
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = global->temperature;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}
	
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }

    //
    // Hierarchical prior
    //
    if (hierarchical_prior != nullptr) {

      if (hierarchical_prior->step() < 0) {
	fprintf(stderr, "error: failed to do hierarchical prior step\n");
	return -1;
      }

      if (chain_rank == 0) {

	chain_history_change_t step;
	
	hierarchical_prior->get_last_step(&step);
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = global->temperature;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}
	
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }      

    //
    // Hierarchical prior
    //
    if (hierarchical_prior != nullptr) {

      if (hierarchical_prior->step() < 0) {
	fprintf(stderr, "error: failed to do hierarchical prior step\n");
	return -1;
      }

      if (chain_rank == 0) {

	chain_history_change_t step;
	
	hierarchical_prior->get_last_step(&step);
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = global->temperature;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}
	
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }      

    //
    // Hierarchical prior
    //
    if (hierarchical_prior != nullptr) {

      if (hierarchical_prior->step() < 0) {
	fprintf(stderr, "error: failed to do hierarchical prior step\n");
	return -1;
      }

      if (chain_rank == 0) {

	chain_history_change_t step;
	
	hierarchical_prior->get_last_step(&step);
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = global->temperature;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}
	
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }      
 
    //
    // Hierarchical prior
    //
    if (hierarchical_prior != nullptr) {

      if (hierarchical_prior->step() < 0) {
	fprintf(stderr, "error: failed to do hierarchical prior step\n");
	return -1;
      }

      if (chain_rank == 0) {

	chain_history_change_t step;
	
	hierarchical_prior->get_last_step(&step);
	
	step.header.likelihood = global->current_likelihood;
	step.header.temperature = global->temperature;
	step.header.hierarchical = global->hierarchical->getparameter(0);
	
	if (chain_history_full(global->ch)) {
	  
	  /*
	   * Flush chain history to file
	   */
	  if (chain_history_write(global->ch,
				  (ch_write_t)fwrite,
				  fp_ch) < 0) {
	    ERROR("error: failed to write chain history segment to file\n");
	    return -1;
	  }
	  
	  if (chain_history_reset(global->ch) < 0) {
	    ERROR("error: failed to reset chain history\n");
	    return -1;
	  }
	  
	}
	
	if (chain_history_add_step(global->ch, &step) < 0) {
	  ERROR("error: failed to add step to chain history\n");
	  return -1;
	}
      }
    }      

   
    if (ptexchange != nullptr && exchange_rate > 0 && ((i + 1) % exchange_rate == 0)) {
	
      int exchanged = ptexchange->step();
      if (exchanged < 0) {
	ERROR("Failed to do PT exchange\n");
	return -1;
      }

      if (exchanged) {
	global->invalidate_residuals();
      }
      
      if (!posteriork && chain_rank == 0 && exchanged == 1) {
	//
	// Flush and reinitialize chain history to deal with completely new model.
	//
	if (chain_history_write(global->ch,
				(ch_write_t)fwrite,
				fp_ch) < 0) {
	  ERROR("error: failed to write chain history segment to file\n");
	  return -1;
	}
	
	if (chain_history_initialise(global->ch,
				     wavetree3d_sub_get_S_v(global->wt),
				     global->current_likelihood,
				     global->temperature,
				     global->hierarchical->getparameter(0)) < 0) {
	  ERROR("error: failed to initialise chain history\n");
	  return -1;
	}
      }

      
    }

    if (chain_rank == 0 && verbosity > 0 && (i + 1) % verbosity == 0) {

      INFO("%03d %6d: %f(%f) %d dc %f lambda %f T %f:",
	   chain_id,
	   i + 1,
	   global->current_likelihood,
	   global->current_log_normalization,
	   current_k,
	   wavetree3d_sub_dc(global->wt),
	   global->hierarchical->getparameter(0),
	   temperature);

      INFO(birth->write_long_stats().c_str());
      INFO(death->write_long_stats().c_str());
      INFO(value->write_long_stats().c_str());

      if (hierarchical != nullptr) {
	INFO(hierarchical->write_long_stats().c_str());
      }

      if (hierarchical_prior != nullptr) {
	INFO(hierarchical_prior->write_long_stats().c_str());
      }
      
      if (ptexchange != nullptr) {
	INFO(ptexchange->write_long_stats().c_str());
      }
    }
  }

  if (chain_rank == 0) {
    std::string filename = mkfilenamerank(output_prefix, "khistogram.txt", chain_id);
    FILE *fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
      ERROR("error: failed to create khistogram file\n");
      return -1;
    }
    for (int i = 0; i < kmax; i ++) {
      fprintf(fp, "%d %d\n", i + 1, khistogram[i]);
    }
    fclose(fp);
    
    if (!posteriork) {
      /*
       * If there are remaining steps to save
       */
      if (chain_history_nsteps(global->ch) > 1) {
	/*
	 * Flush chain history to file
	 */
	if (chain_history_write(global->ch,
				(ch_write_t)fwrite,
				fp_ch) < 0) {
	  ERROR("error: failed to write chain history segment to file\n");
	  return -1;
	}
      }
      fclose(fp_ch);
    }
    
    filename = mkfilenamerank(output_prefix, "acceptance.txt", chain_id);
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
      ERROR("error: failed to create acceptance file\n");
      return -1;
    }
    fprintf(fp, birth->write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, death->write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, value->write_long_stats().c_str());
    fprintf(fp, "\n");
    if (hierarchical != nullptr) {
      fprintf(fp, hierarchical->write_long_stats().c_str());
      fprintf(fp, "\n");
    }
    if (ptexchange != nullptr) {
      fprintf(fp, ptexchange->write_long_stats().c_str());
      fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    filename = mkfilenamerank(output_prefix, "final_model.txt", chain_id);
    if (wavetree3d_sub_save(global->wt, filename.c_str()) < 0) {
      ERROR("error: failed to save final model\n");
      return -1;
    }

    filename = mkfilenamerank(output_prefix, "residuals.txt", chain_id);
    int nres = global->get_residual_size();
    const double *res = global->get_mean_residuals();
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
      ERROR("Failed to create residuals file");
      return -1;
    }
    for (int i = 0; i < nres; i ++) {
      fprintf(fp, "%.9g\n", res[i]);
    }
    fclose(fp);

    filename = mkfilenamerank(output_prefix, "residuals_normed.txt", chain_id);
    res = global->get_mean_normed_residuals();
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
      ERROR("Failed to create residuals squared file");
      return -1;
    }
    for (int i = 0; i < nres; i ++) {
      fprintf(fp, "%.9g\n", res[i]);
    }
    fclose(fp);

    filename = mkfilenamerank(output_prefix, "residuals_hist.txt", chain_id);
    if (!global->save_residual_histogram(filename.c_str())) {
      ERROR("Failed to save residual histogram");
      return -1;
    }
    
    filename = mkfilenamerank(output_prefix, "residuals_cov.txt", chain_id);
    if (!global->save_residual_covariance(filename.c_str())) {
      ERROR("Failed to save residual covariance");
      return -1;
    }

    filename = mkfilenamerank(output_prefix, "zoffset.txt", chain_id);
    if (!save_zoffset(filename.c_str(), global->zoffset, global->depth)) {
      ERROR("Failed to save zoffset");
      return -1;
    }
    
  }
  
  MPI_Finalize();
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <file>               Input observations file\n"
	  " -I|--initial <file>             Starting model file\n"
	  " -M|--prior-file <file>          Prior/Proposal file\n"
	  " -H|--hierarchical <filename>    Hierarchical prior/proposal file\n"
	  " -L|--hierarchical-prior-sigma <float> Hierarchical prior proposal std dev\n"
	  " -o|--output <path>              Output prefix for output files\n"
	  "\n"
	  " -x|--degree-x <int>             Number of samples in x/lon direction as power of 2\n"
	  " -y|--degree-y <int>             Number of samples in y/lat direction as power of 2\n"
	  " -z|--degree-z <int>             Number of frequency samples as power of 2\n"
	  "\n"
	  " -t|--total <int>                Total number of iterations\n"
	  " -S|--seed <int>                 Random number seed\n"
	  "\n"
	  " -l|--initial-lambda <float>         Std deviation for lambda scaling sampling\n"
	  "\n"
	  " -k|--kmax <int>                 Max. no. of coefficients\n"
	  "\n"
	  " -B|--birth-probability <float>  Birth probability\n"
	  " -P|--posteriork                 Posterior k simulation\n"
	  "\n"
	  " -w|--wavelet-xy <int>           Wavelet basis to use for lateral direction\n"
	  " -W|--wavelet-z <int>            Wavelet basis to use for depth/frequency direction\n"
	  "\n"
	  " -c|--chains <int>               No. of chains per temperature\n"
	  " -T|--temperatures <int>         No. of temperature levels\n"
	  " -m|--max-temperature <float>    Max. Temperature\n"
	  " -e|--exchange-rate <int>        No. of steps between exchange proposals\n"
	  "\n"
	  " -v|--verbosity <int>            Number steps between status printouts (0 = disable)\n"
	  " -h|--help                       Show usage information\n"
	  "\n",
	  pname);
}
