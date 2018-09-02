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

#include <string.h>

#include "wavetomo2dexception.hpp"
#include "wavetomo2dobservations.hpp"
#include "wavetomo2dutil.hpp"

#include "globalslice.hpp"

extern "C" {
  #include "hnk_cartesian_nonsquare.h"

  #include "slog.h"
};

static constexpr double LARGE_LIKELIHOOD = 1e99;
static constexpr int SUBTILE = 1;

const int CHAIN_STEPS = 1000000;

int globalslice_coordtoindex(void *user, int i, int j, int k, int depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  return wavetree2d_sub_from_2dindices(wt, i, j);
}

int globalslice_indextocoord(void *user, int index, int *i, int *j, int *k, int *depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  if (wavetree2d_sub_2dindices(wt, index, i, j) < 0) {
    return -1;
  }

  k = 0;
  *depth = wavetree2d_sub_depthofindex(wt, index);
  return 0;
}

GlobalSlice::GlobalSlice(const char *filename,
			 const char *initial_model,
			 const char *prior_file,
			 int _degreex,
			 int _degreey,
			 int _degreez,
			 int super_resolution,
			 int _slice,
			 double _lonmin,
			 double _lonmax,
			 double _latmin,
			 double _latmax,
			 int seed,
			 int _kmax,
			 double lambda,
			 bool _posteriork,
			 int xywavelet,
			 bool _linear) :
  kmax(_kmax),
  treemaxdepth(-1),
  linear(_linear),
  wt(nullptr),
  ch(nullptr),
  hnk(nullptr),
  proposal(nullptr),
  degreex(_degreex),
  degreey(_degreey),
  degreez(_degreez),
  slice(_slice),
  lonmin(_lonmin),
  lonmax(_lonmax),
  latmin(_latmin),
  latmax(_latmax),
  observations(nullptr),
  model(nullptr),
  workspace(nullptr),
  mean_residual_n(0),
  residual(nullptr),
  mean_residual(nullptr),
  last_valid_residual(nullptr),
  residual_normed(nullptr),
  mean_residual_normed(nullptr),
  last_valid_residual_normed(nullptr),
  residuals_valid(false),
  residual_hist_bins(100),
  residual_hist_min(-5.0),
  residual_hist_max(5.0),
  residual_hist(nullptr),
  width(-1),
  height(-1),
  size(-1),
  ncoeff(-1),
  zoffset(nullptr),
  hierarchical(nullptr),
  current_likelihood(-1.0),
  coeff_hist(nullptr),
  random(seed),
  posteriork(_posteriork),
  xywaveletf(nullptr),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1),
  temperature(1.0),
  trace_offsets(nullptr),
  trace_sizes(nullptr),
  residual_offsets(nullptr),
  residual_sizes(nullptr),
  cov_n(-1),
  super_image_width(-1),
  super_image_height(-1),
  super_image_size(-1),
  dLdI(nullptr),
  wtmapper(nullptr)
{
  if (degreex < 0 || degreex >= 16 ||
      degreey < 0 || degreey >= 16 ||
      degreez < 0 || degreez >= 16) {
    throw WAVETOMO2DEXCEPTION("Degree(s) out of range: %d x %d x %d\n", degreex, degreey, degreez);
  }
  
  xywaveletf = wavelet_inverse_function_from_id(xywavelet);
  if (xywaveletf == nullptr) {
    throw WAVETOMO2DEXCEPTION("Invalid horizontal wavelet %d\n", xywavelet);
  }

  if (!posteriork) {
    //
    // Load observations
    //
    observations = new wavetomo2dobservations<LonLat<>>(degreex, degreey, degreez,
							super_resolution,
							lonmin, lonmax,
							latmin, latmax,
							filename,
							false, // Travel time
							linear,
							xywaveletf);

    zoffset = new double[observations->get_frequency_count()];

    observations->compute_mean_velocities(zoffset, observations->get_frequency_count());

    INFO("%d stations, %d traces, %d frequencies, %d observations\n",
	 observations->get_station_count(),
	 observations->get_trace_count(),
	 observations->get_frequency_count(),
	 observations->get_observation_count());

    for (int i = 0; i < (int)observations->get_frequency_count(); i ++) {
      INFO("%2d %10.6f\n", i, zoffset[i]);
    }

  }

  hierarchical = new independentgaussianhierarchicalmodel();
  hierarchical->setparameter(0, lambda);
  
  wt = wavetree2d_sub_create(degreex, degreey, 0.0);
  if (wt == NULL) {
    throw WAVETOMO2DEXCEPTION("Failed to create wavetree\n");
  }

  width = wavetree2d_sub_get_width(wt);
  height = wavetree2d_sub_get_height(wt);
  slicestride = width * height;
  size = wavetree2d_sub_get_size(wt);
  ncoeff = wavetree2d_sub_get_ncoeff(wt);
  treemaxdepth = wavetree2d_sub_maxdepth(wt);

  INFO("Image: %d x %d\n", width, height);

  if (!posteriork) {

    model = new double[size];
    int workspacesize = width;
    if (height > workspacesize) {
      workspacesize = height;
    }
    workspace = new double[workspacesize];

    int ntotal = observations->get_slice_observation_count(slice);
    INFO("Data: %d total points\n", ntotal);

    residual_size = ntotal;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double [residual_size];

    residual_hist = new int[residual_size * residual_hist_bins];

    reset_residuals();
  } else {
    residual_size = 0;
    residual = nullptr;
    mean_residual = nullptr;
    last_valid_residual = nullptr;
    residual_normed = nullptr;
    mean_residual_normed = nullptr;
    last_valid_residual_normed = nullptr;
    residual_hist = nullptr;
  }
  
  if (initial_model == NULL) {
    if (wavetree2d_sub_initialize(wt, 0.0) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to initialize wavetree\n");
    }
  } else {

    if (wavetree2d_sub_load(wt, initial_model) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to load initial model: %s\n", initial_model);
    }

    INFO("Loaded model with %d coefficients\n", wavetree2d_sub_coeff_count(wt));
  }

  //
  // Hnk Ratio
  //
  if (kmax > ncoeff) {
    INFO("Warning: kmax truncated to %d\n", ncoeff);
    kmax = ncoeff;
  }
  
  hnk = hnk_cartesian_nonsquare_2D_create_sub(degreex,
					      degreey,
					      kmax);
  if (hnk == NULL) {
    throw WAVETOMO2DEXCEPTION("Failed to create hnk table\n");
  }

  //
  // Chain History
  //
  ch = chain_history_create(CHAIN_STEPS);
  if (ch == nullptr) {
    throw WAVETOMO2DEXCEPTION("Failed to create chain history\n");
  }

  
  //
  // Initialse coeff histogram
  //
  coeff_hist = coefficient_histogram_create(ncoeff, 100, -1.0, 1.0,
					    globalslice_coordtoindex,
					    globalslice_indextocoord,
					    wt);
  if (coeff_hist == NULL) {
    throw WAVETOMO2DEXCEPTION("Failed to create coefficient histogram\n");
  }

  //
  // Create proposal structure.
  //
  if (prior_file != nullptr) {
    proposal = wavetree_pp_load(prior_file, seed, coeff_hist);
    if (proposal == NULL) {
      throw WAVETOMO2DEXCEPTION("Failed to load proposal file\n");
    }
    
    for (int i = 0; i < ncoeff; i ++) {
      int coeffdepth = wavetree2d_sub_depthofindex(wt, i);
      
      int ii, ij;
      if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to get 2d indices\n");
      }
      
      double vmin, vmax;
      if (wavetree_pp_prior_range2d(proposal,
				    ii,
				    ij,
				    coeffdepth,
				    treemaxdepth,
				    0.0,
				    &vmin,
				    &vmax) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to get coefficient range\n");
      }
      
      if (coefficient_histogram_set_range(coeff_hist, 
					  i,
					  vmin,
					  vmax) < 0) {
	throw WAVETOMO2DEXCEPTION("Failed to set coefficient histogram range\n");
      }
    }
  }
    
}

GlobalSlice::~GlobalSlice()
{
  delete hierarchical;
  delete observations;
  delete [] zoffset;

  delete [] model;
  delete [] workspace;

  delete [] residual;
  delete [] mean_residual;
  delete [] last_valid_residual;
  delete [] residual_normed;
  delete [] mean_residual_normed;
  delete [] last_valid_residual_normed;
  delete [] residual_hist;

  coefficient_histogram_destroy(coeff_hist);
  hnk_destroy(hnk);
  wavetree2d_sub_destroy(wt);
  wavetree_pp_destroy(proposal);
}

double
GlobalSlice::image_likelihood(const double *image_model,
			      double &log_normalization)
{
    log_normalization = 0.0;
    return observations->single_frequency_likelihood(slice,
						     image_model,
						     0.0,
						     hierarchical,
						     residual,
						     residual_normed,
						     log_normalization);
}
			 
double
GlobalSlice::likelihood(double &log_normalization)
{
  if (!posteriork) {
    
    //
    // Get tree model wavelet coefficients
    //
    memset(model, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, model, size) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(model,
    			       width,
    			       height,
    			       width,
    			       workspace,
    			       xywaveletf,
    			       xywaveletf,
    			       SUBTILE) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    log_normalization = 0.0;
    return observations->single_frequency_likelihood(slice,
						     model,
						     zoffset[slice],
						     hierarchical,
						     residual,
						     residual_normed,
						     log_normalization);
    
  } else {
    return 1.0;
  }
}

double
GlobalSlice::likelihood_gradient(double &log_normalization)
{
  //
  // Get tree model wavelet coefficients
  //
  memset(model, 0, sizeof(double) * size);
  if (wavetree2d_sub_map_to_array(wt, model, size) < 0) {
    throw WAVETOMO2DEXCEPTION("Failed to map model to array\n");
  }
  
  //
  // Inverse wavelet transform
  //
  if (generic_lift_inverse2d(model,
			     width,
			     height,
			     width,
			     workspace,
			     xywaveletf,
			     xywaveletf,
			     SUBTILE) < 0) {
    throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
  }
  
  log_normalization = 0.0;


  // Normal likelihood calculation but with dLdI (gradient of likelihood wrt model image)
  memset(dLdI, 0, sizeof(double) * super_image_size);
  double like = observations->single_frequency_likelihood_gradient(slice,
								   model,
								   dLdI,
								   zoffset[slice],
								   hierarchical,
								   residual,
								   residual_normed,
								   log_normalization);

  //
  // Back project dLdI to TDT wavelet coefficients
  //
  wtmapper->backproject(wt,
			width,
			height,
			xywaveletf,
			xywaveletf,
			workspace,
			observations,
			zoffset[slice],
			dLdI);

  return like;
}



double
GlobalSlice::hierarchical_likelihood(double proposed_lambda_scale,
				     double &log_normalization)
{
  log_normalization = 0.0;

  if (!posteriork) {

    if (!residuals_valid) {
      double x;
      (void)likelihood(x);
      accept();
    }

    log_normalization = 0.0;

    double oldlambda = hierarchical->getparameter(0);
    hierarchical->setparameter(0, proposed_lambda_scale);
    double sum = observations->single_frequency_hierarchical_likelihood(slice,
									hierarchical,
									last_valid_residual,
									last_valid_residual_normed,
									log_normalization);
    hierarchical->setparameter(0, oldlambda);

    return sum;
    
  } else {
    return 1.0;
  }
}

void
GlobalSlice::initialize_mpi(MPI_Comm _communicator, double _temperature)
{
  MPI_Comm_dup(_communicator, &communicator);

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }
  
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw WAVETOMO2DEXCEPTION("MPI Failure\n");
  }

  INFO("GlobalSlice MPI: %3d/%3d", mpi_rank, mpi_size);
  
  trace_offsets = new int[mpi_size];
  trace_sizes = new int[mpi_size];
  residual_offsets = new int[mpi_size];
  residual_sizes = new int[mpi_size];

  int traces = 0;
  if (observations != nullptr) {
    traces = observations->get_trace_count();
  }
  
  int processes = mpi_size;

  //
  // Evenly distribute traces amoungst processes
  //
  for (int i = 0; i < mpi_size; i ++) {
    trace_sizes[i] = traces/processes;

    traces -= trace_sizes[i];
    processes --;
  }

  trace_offsets[0] = 0;
  residual_offsets[0] = 0;

  residual_sizes[0] = 0;
  for (int i = 0; i < trace_sizes[0]; i ++) {
    if (observations != nullptr) {
      residual_sizes[0] += observations->get_trace_slice_observation_count(i);
    }
  }

  INFO("Split: %4d %4d %4d %4d", trace_offsets[0], trace_sizes[0], residual_offsets[0], residual_sizes[0]);
  
  for (int i = 1; i < mpi_size; i ++) {
    trace_offsets[i] = trace_offsets[i - 1] + trace_sizes[i - 1];
    residual_offsets[i] = residual_offsets[i - 1] + residual_sizes[i - 1];

    residual_sizes[i] = 0;
    for (int j = 0; j < trace_sizes[i]; j ++) {
      if (observations != nullptr) {
	residual_sizes[i] += observations->get_trace_slice_observation_count(trace_offsets[i] + j);
      }
    }
								     
    INFO("Split: %4d %4d %4d %4d", trace_offsets[i], trace_sizes[i], residual_offsets[i], residual_sizes[i]);
  }

  if (observations != nullptr) {
    if (trace_offsets[mpi_size - 1] + trace_sizes[mpi_size - 1] != observations->get_trace_count()) {
      throw WAVETOMO2DEXCEPTION("Trace sharing intialization failure (trace mismatch %d + %d != %d)",
				trace_offsets[mpi_size - 1],
				trace_sizes[mpi_size - 1],
				observations->get_trace_count());
    }
    
    if (residual_offsets[mpi_size - 1] + residual_sizes[mpi_size - 1] != observations->get_slice_observation_count(slice)) {
      throw WAVETOMO2DEXCEPTION("Trace sharing intialization failure (residual mismatch %d + %d != %d)",
				residual_offsets[mpi_size - 1],
				residual_sizes[mpi_size - 1],
				observations->get_observation_count());
    }
  }
  
  temperature = _temperature;
}

double
GlobalSlice::likelihood_mpi(double &log_normalization)
{
  if (communicator == MPI_COMM_NULL || mpi_rank < 0 || mpi_size < 0) {
    throw WAVETOMO2DEXCEPTION("MPI Parameters unset\n");
  }
  
  if (!posteriork) {
    
    //
    // Get tree model wavelet coefficients
    //
    memset(model, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, model, size) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(model,
			       width,
			       height,
			       width,
			       workspace,
			       xywaveletf,
			       xywaveletf,
			       SUBTILE) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }
    
    double local_sum = 0.0;
    double local_log_normalization = 0.0;
    
    double *res = residual + residual_offsets[mpi_rank];
    double *resn = residual_normed + residual_offsets[mpi_rank];
    
    for (int mi = 0, i = trace_offsets[mpi_rank]; mi < trace_sizes[mpi_rank]; mi ++, i ++) {
      
      local_sum += observations->single_frequency_trace_likelihood(slice,
								   i,
								   model,
								   zoffset[slice],
								   hierarchical,
								   res,
								   resn,
								   local_log_normalization);
      
      res += observations->get_trace_slice_observation_count(i);
      resn += observations->get_trace_slice_observation_count(i);
    }
    
    double total;

    if (MPI_Reduce(&local_log_normalization, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Likelihood failed in broadcast\n");
    }
    
    log_normalization = total;
    
    if (MPI_Reduce(&local_sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw WAVETOMO2DEXCEPTION("Likelihood failed in broadcast\n");
    }

    
    if (mpi_rank == 0) {
      MPI_Gatherv(MPI_IN_PLACE,
		  residual_sizes[mpi_rank],
		  MPI_DOUBLE,
		  residual,
		  residual_sizes,
		  residual_offsets,
		  MPI_DOUBLE,
		  0,
		  communicator);
    } else {
      MPI_Gatherv(residual + residual_offsets[mpi_rank],
		  residual_sizes[mpi_rank],
		  MPI_DOUBLE,
		  residual,
		  residual_sizes,
		  residual_offsets,
		  MPI_DOUBLE,
		  0,
		  communicator);
    }

    if (mpi_rank == 0) {
      MPI_Gatherv(MPI_IN_PLACE,
		  residual_sizes[mpi_rank],
		  MPI_DOUBLE,
		  residual_normed,
		  residual_sizes,
		  residual_offsets,
		  MPI_DOUBLE,
		  0,
		  communicator);
    } else {
      MPI_Gatherv(residual_normed + residual_offsets[mpi_rank],
		  residual_sizes[mpi_rank],
		  MPI_DOUBLE,
		  residual_normed,
		  residual_sizes,
		  residual_offsets,
		  MPI_DOUBLE,
		  0,
		  communicator);
    }
    
    return total;
    
  } else {
    log_normalization = 0.0;
    
    return 1.0;
  }
}

double
GlobalSlice::hierarchical_likelihood_mpi(double proposed_lambda_scale,
					 double &log_normalization)
{
  double like;

  log_normalization = 0.0;

  if (!residuals_valid) {
    double x;
    like = likelihood_mpi(x);
    accept();
  }
  
  //
  // Since this is just a sum over residuals, just compute on root node and broadcast.
  //
  if (mpi_rank == 0) {
    like = hierarchical_likelihood(proposed_lambda_scale, log_normalization);
  }
  
  MPI_Bcast(&like, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&log_normalization, 1, MPI_DOUBLE, 0, communicator);

  return like;
}

void
GlobalSlice::reset_residuals()
{
  mean_residual_n = 0;
  for (int i = 0; i < residual_size; i ++) {
    residual[i] = 0.0;
    mean_residual[i] = 0.0;
    last_valid_residual[i] = 0.0;
      
    residual_normed[i] = 0.0;
    mean_residual_normed[i] = 0.0;
    last_valid_residual_normed[i] = 0.0;

    for (int j = 0; j < residual_hist_bins; j ++) {
      residual_hist[i * residual_hist_bins + j] = 0;
    }
  }

  cov_n = 0;

  for (int i = 0; i < (int)cov_count.size(); i ++) {
    int N = cov_count[i];

    for (int j = 0; j < N; j ++) {
      cov_delta[i][j] = 0.0;
      cov_mu[i][j] = 0.0;
    }

    N = N*N;
    for (int j = 0; j < N; j ++) {
      cov_sigma[i][j] = 0.0;
    }
  }
}

void
GlobalSlice::invalidate_residuals()
{
  residuals_valid = false;
}

void
GlobalSlice::accept()
{
  residuals_valid = true;
  if (!posteriork) {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residual[i] = residual[i];
      last_valid_residual_normed[i] = residual_normed[i];
    }

    update_residual_mean();
    update_residual_covariance();
  }
}

void
GlobalSlice::accept_hierarchical()
{
}

void
GlobalSlice::reject()
{
  update_residual_mean();
}

void
GlobalSlice::reject_hierarchical()
{
}

void
GlobalSlice::update_residual_mean()
{
  mean_residual_n ++;

  for (int i = 0; i < residual_size; i ++) {

    double delta = last_valid_residual[i] - mean_residual[i];
    mean_residual[i] += delta/(double)(mean_residual_n);

    delta = last_valid_residual_normed[i] - mean_residual_normed[i];
    mean_residual_normed[i] += delta/(double)(mean_residual_n);

    int hi = (int)((last_valid_residual_normed[i] - residual_hist_min)/(residual_hist_max - residual_hist_min) * (double)residual_hist_bins);
    if (hi >= 0 && hi < residual_hist_bins) {
      residual_hist[i * residual_hist_bins + hi] ++;
    }
  }
}

void
GlobalSlice::update_residual_covariance()
{
  double *p = last_valid_residual;

  for (int k = 0; k < residual_size; k ++) {

    cov_n ++;
    
    for (int i = 0; i < (int)cov_count.size(); i ++) {

      int N = cov_count[i];
      
      for (int j = 0; j < N; j ++) {
  	cov_delta[i][j] = (p[j] - cov_mu[i][j])/(double)(cov_n);
  	cov_mu[i][j] += cov_delta[i][j];
      }

      for (int j = 0; j < N; j ++) {
  	for (int l = j; l < N; l ++) {

  	  cov_sigma[i][j * N + l] +=
  	    (double)(cov_n - 1)*cov_delta[i][j]*cov_delta[i][l] -
  	    cov_sigma[i][j * N + l]/(double)(cov_n);
  	}
      }

      p += N;
    }
  }
}


int
GlobalSlice::get_residual_size() const
{
  return residual_size;
}

const double *
GlobalSlice::get_mean_residuals() const
{
  return mean_residual;
}
  
const double *
GlobalSlice::get_mean_normed_residuals() const
{
  return mean_residual_normed;
}

bool
GlobalSlice::save_residuals(const char *filename)
{
  if (observations != nullptr) {
    return observations->save_residuals(slice, filename, mean_residual, mean_residual_normed);
  } else {
    return true;
  }
}

bool
GlobalSlice::save_residual_histogram(const char *filename) const
{
  if (observations != nullptr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      ERROR("Failed to create file");
      return false;
    }
    
    fprintf(fp, "%d %d %f %f\n", residual_size, residual_hist_bins, residual_hist_min, residual_hist_max);
    for (int i = 0; i < residual_size; i ++) {
      for (int j = 0; j < residual_hist_bins; j ++) {
	fprintf(fp, "%d ", residual_hist[i * residual_hist_bins + j]);
      }
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
  
  return true;
}

bool
GlobalSlice::save_residual_covariance(const char *filename) const
{
  if (observations != nullptr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      ERROR("Failed to create file\n");
      return false;
    }
    
    fprintf(fp, "%d\n", (int)cov_count.size());
    
    for (int i = 0; i < (int)cov_count.size(); i ++) {
      
      int N = cov_count[i];
      
      fprintf(fp, "%d\n", N);
      for (int j = 0; j < N; j ++) {
	fprintf(fp, "%.9g ", cov_mu[i][j]);
      }
      fprintf(fp, "\n");
      
      for (int j = 0; j < N; j ++) {
	for (int k = 0; k < N; k ++) {
	  fprintf(fp, "%.9g ", cov_sigma[i][j * N + k]);
	}
	fprintf(fp, "\n");
      }
    }
    
    fclose(fp);
  }

  return true;
}

bool
GlobalSlice::save_zoffset(const char *filename, int slice)
{
  if (observations != nullptr) {
    return ::save_zoffset(filename, zoffset + slice, 1);
  } else {
    return true;
  }
}

generic_lift_inverse1d_step_t
GlobalSlice::wavelet_inverse_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_inverse1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_inverse1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_inverse1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_inverse1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_inverse1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_inverse1d_cdf97_step;

  default:
    return nullptr;
  }
}

generic_lift_forward1d_step_t
GlobalSlice::wavelet_forward_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_forward1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_forward1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_forward1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_forward1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_forward1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_forward1d_cdf97_step;

  default:
    return nullptr;
  }
}

void GlobalSlice::set_max_depth(int md)
{
  int maxd = wavetree2d_sub_maxdepth(wt);
  if (md > maxd) {
    md = maxd;
  }

  if (md < 1) {
    md = 1;
  }
  
  treemaxdepth = md;

  int dmaxk = hnk_get_maxk_at_h(hnk, treemaxdepth);

  printf("Depth maxk: %d\n", dmaxk);
  if (kmax > dmaxk) {
    kmax = dmaxk;
  }
  printf("New   maxk: %d\n", kmax);
  

}

bool
GlobalSlice::initialize_gradient()
{
  super_image_size = observations->image_size(super_image_width, super_image_height);

  dLdI = new double[super_image_size];
  wtmapper = new WavetreeMapper(kmax + 1);
  
  return true;
}
