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
#ifndef wavetomo2dobservations_hpp
#define wavetomo2dobservations_hpp

#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>

#include "traveltimefield.hpp"
#include "velocityfield.hpp"

#include "linearweights.hpp"

#include "wavetomo2dexception.hpp"
#include "rng.hpp"
#include "hierarchicalmodel.hpp"

extern "C" {
#include "generic_lift.h"
};

template <
  typename coordinate
>
class wavetomo2dobservations {
public:

  static constexpr double RADIUS = 6371.0;
  static constexpr double LARGE_LIKELIHOOD = 1e99;

  //
  // Constructor for synthetic
  //
  wavetomo2dobservations(int degreex,
			 int degreey,
			 int degreez,
			 int _super_resolution,
			 double _xmin,
			 double _xmax,
			 double _ymin,
			 double _ymax,
			 bool _traveltime,
			 bool _linear,
			 generic_lift_inverse1d_step_t _winverse) :
    width(1 << degreex),
    height(1 << degreey),
    super_resolution(_super_resolution),
    super_width(-1),
    super_height(-1),
    super_size(-1),
    super_image(nullptr),
    super_workspace(nullptr),
    xmin(_xmin),
    xmax(_xmax),
    ymin(_ymin),
    ymax(_ymax),
    velocity(nullptr),
    refinement(1),
    traveltime(_traveltime),
    linear(_linear),
    winverse(_winverse)
  {
    if (super_resolution < 0) {
      throw WAVETOMO2DEXCEPTION("Invalid super resolution %d\n", super_resolution);
    }

    if (super_resolution > 0) {
      super_width = width << super_resolution;
      super_height = height << super_resolution;
    } else {
      super_width = width;
      super_height = height;
    }

    super_size = super_width * super_height;
      
    int super_workspace_size = super_width;
    if (super_height > super_width) {
      super_workspace_size = super_height;
    }
      
    super_image = new double[super_width * super_height];
    super_workspace = new double[super_workspace_size];

    velocity = new VelocityField<double>(super_width, super_height);
  }

  //
  // Constructor for inversion
  //
  wavetomo2dobservations(int degreex,
			 int degreey,
			 int degreez,
			 int _super_resolution,
			 double _xmin,
			 double _xmax,
			 double _ymin,
			 double _ymax,
			 const char *filename,
			 bool _traveltime,
			 bool _linear,
			 generic_lift_inverse1d_step_t _winverse) :
    width(1 << degreex),
    height(1 << degreey),
    super_resolution(_super_resolution),
    super_width(-1),
    super_height(-1),
    super_size(-1),
    super_image(nullptr),
    super_workspace(nullptr),
    xmin(_xmin),
    xmax(_xmax),
    ymin(_ymin),
    ymax(_ymax),
    velocity(nullptr),
    refinement(1),
    traveltime(_traveltime),
    linear(_linear),
    winverse(_winverse)
  {
    if (super_resolution < 0) {
      throw WAVETOMO2DEXCEPTION("Invalid super resolution %d\n", super_resolution);
    }

    if (super_resolution > 0) {
      super_width = width << super_resolution;
      super_height = height << super_resolution;
    } else {
      super_width = width;
      super_height = height;
    }

    super_size = super_width * super_height;
      
    int super_workspace_size = super_width;
    if (super_height > super_width) {
      super_workspace_size = super_height;
    }
      
    super_image = new double[super_width * super_height];
    super_workspace = new double[super_workspace_size];

    velocity = new VelocityField<double>(super_width, super_height);
    
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
      throw WAVETOMO2DEXCEPTION("Failed to open file\n");
    }
    
    //
    // Read the stations
    //
    int nstations;
    if (fscanf(fp, "%d\n", &nstations) != 1) {
      throw WAVETOMO2DEXCEPTION("Failed to read no. stations");
    }

    station_xmin = 1e9;
    station_xmax = -1e9;
    station_ymin = 1e9;
    station_ymax = -1e9;
    
    for (int i = 0; i < nstations; i ++) {
      char code[256];
      double x, y;

      if (fscanf(fp, "%s %lf %lf\n", code, &x, &y) != 3) {
	throw WAVETOMO2DEXCEPTION("Failed to read station");
      }

      if (x < xmin || x > xmax) {
	throw WAVETOMO2DEXCEPTION("Station out of range x: %s %f (%f %f)\n", code, x, xmin, xmax);
      }

      if (y < ymin || y > ymax) {
	throw WAVETOMO2DEXCEPTION("Station out of range y: %s %f (%f %f)\n", code, y, ymin, ymax);
      }

      if (x < station_xmin) {
	station_xmin = x;
      }
      if (x > station_xmax) {
	station_xmax = x;
      }
      if (y < station_ymin) {
	station_ymin = y;
      }
      if (y > station_ymax) {
	station_ymax = y;
      }
      
      stations.push_back(station(code, coordinate(x, y)));

    }

    //
    // Read the frequencies
    //
    int nfrequencies;
    if (fscanf(fp, "%d\n", &nfrequencies) != 1) {
      throw WAVETOMO2DEXCEPTION("Failed to read no. frequencies");
    }

    for (int i = 0; i < nfrequencies; i ++) {
      double freq;

      if (fscanf(fp, "%lf\n", &freq) != 1) {
	throw WAVETOMO2DEXCEPTION("Failed to read frequency");
      }

      frequencies.push_back(freq);

    }
    
    //
    // Read the optimized tracing and build the empty trace structures.
    //
    int ntraces;
    if (fscanf(fp, "%d\n", &ntraces) != 1) {
      throw WAVETOMO2DEXCEPTION("Failed to read no. traces");
    }

    for (int i = 0; i < ntraces; i ++) {
      char codeA[256];
      int npairs;
      coordinate cA;
      
      if (fscanf(fp, "%s %d", codeA, &npairs) != 2) {
	throw WAVETOMO2DEXCEPTION("Failed to read trace header");
      }

      if (!find_station_coordinate(codeA, cA)) {
	throw WAVETOMO2DEXCEPTION("Failed to find coordinate for station %s", codeA);
      }

      trace *t = new trace(codeA, cA, velocity, xmin, xmax, ymin, ymax, refinement);
      traces.push_back(t);
	
      for (int j = 0; j < npairs; j ++) {

	char codeB[256];
	if (fscanf(fp, "%s", codeB) != 1) {
	  throw WAVETOMO2DEXCEPTION("Failed to read trace pair");
	}

	coordinate cB;
	if (!find_station_coordinate(codeB, cB)) {
	  throw WAVETOMO2DEXCEPTION("Failed to find coordinate for station %s", codeB);
	}

	(void)t->add_tracepair(codeB,
			       cB,
			       nfrequencies);
      }
    }
    
    //
    // Read the paths and fill the trace structures
    //
    int npaths;
    if (fscanf(fp, "%d\n", &npaths) != 1) {
      throw WAVETOMO2DEXCEPTION("Failed to read no. paths");
    }

    for (int i = 0; i < npaths; i ++) {

      char codeA[256];
      char codeB[256];
      double distkm;

      if (fscanf(fp, "%s %s %lf\n", codeA, codeB, &distkm) != 3) {
	throw WAVETOMO2DEXCEPTION("Failed to read path header");
      }

      tracepair *tp = find_tracepair(codeA, codeB);
      if (tp == nullptr) {
	throw WAVETOMO2DEXCEPTION("Failed to find tracepair %s - %s", codeA, codeB);
      }

      tp->set_distance(distkm);

      // Read dispersion
      for (int j = 0; j < nfrequencies; j ++) {
	double mean, median, mode, stddev;

	if (fscanf(fp, "%lf %lf %lf %lf\n", &mean, &median, &mode, &stddev) != 4) {
	  throw WAVETOMO2DEXCEPTION("Failed to read dispersion entry");
	}

	//tp->set_frequency_observation(j, mean, median, mode, stddev);
	tp->set_frequency_observation(j, mean, median, mode, stddev);
      }

    }

    fclose(fp);

    if (linear) {
      if (!compute_linear_weights()) {
	throw WAVETOMO2DEXCEPTION("Failed to compute linear weights");
      }
    }
  }
  
  ~wavetomo2dobservations()
  {
    delete velocity;
    
    for (auto t: traces) {
      delete t;
    }

    delete [] super_image;
    delete [] super_workspace;
  }

  bool save(const char *filename)
  {
    FILE *fp = fopen(filename, "w");

    if (fp == NULL) {
      return false;
    }
    
    //
    // Stations
    //
    fprintf(fp, "%d\n", (int)stations.size());

    for (auto &s: stations) {
      fprintf(fp, "%s %15.9f %15.9f\n", s.code.c_str(), s.location.X(), s.location.Y());
    }

    //
    // Frequencies
    //
    fprintf(fp, "%d\n", (int)frequencies.size());

    for (auto &f: frequencies) {
      fprintf(fp, "%15.9f\n", f);
    }

    //
    // Trace schedule
    //
    fprintf(fp, "%d\n", (int)traces.size());
    size_t total_paths = 0;

    for (auto t: traces) {
      fprintf(fp, "%s %d", t->code.c_str(), (int)t->paths.size());
      total_paths += t->paths.size();

      for (auto tp: t->paths) {
	fprintf(fp, " %s", tp->code.c_str());
      }

      fprintf(fp, "\n");
    }

    //
    // Path observations
    //
    fprintf(fp, "%d\n", (int)total_paths);

    for (auto t: traces) {

      for (auto tp: t->paths) {

	fprintf(fp, "%s %s %15.9f\n", t->code.c_str(), tp->code.c_str(), tp->distkm);

	for (size_t fi = 0; fi < tp->nfreq; fi ++) {
	  fprintf(fp, "%15.9f %15.9f %15.9f %15.9f\n", tp->mean[fi], tp->median[fi], tp->mode[fi], tp->stddev[fi]);
	}
      }
    }

    fclose(fp);
    return true;
  }

  int get_station_count() const
  {
    return stations.size();
  }

  int get_frequency_count() const
  {
    return frequencies.size();
  }

  double get_frequency(int i) const
  {
    return frequencies[i];
  }

  int get_trace_count() const
  {
    return traces.size();
  }

  int get_trace_slice_observation_count(size_t trace_index) const
  {
    trace *t = traces[trace_index];

    return t->paths.size();
  }

  int get_trace_observation_count(size_t trace_index) const
  {
    trace *t = traces[trace_index];

    return t->paths.size() * get_frequency_count();
  }
  
  int get_observation_count() const
  {
    int c = 0;
    for (auto t: traces) {
      c += t->paths.size();
    }

    return c * get_frequency_count();
  }

  int get_slice_observation_count(size_t frequency_index) const
  {
    int c = 0;
    for (auto t: traces) {
      c += t->paths.size();
    }

    return c;
  }
  
  void compute_mean_velocities(double *zoffset, int N) const
  {
    int meann;
    double mean;
    double delta;

    if (N != (int)frequencies.size()) {
      throw WAVETOMO2DEXCEPTION("Size mismatch in zoffset: %d != %d", N, (int)frequencies.size());
    }
	
    for (size_t fi = 0; fi < frequencies.size(); fi ++) {
      
      mean = 0.0;
      meann = 0;
      
      for (auto t: traces) {
	
	for (auto tp: t->paths) {
	  
	  if (tp->distkm < 0.0) {
	    throw WAVETOMO2DEXCEPTION("Unset observation");
	  }

	  meann ++;
	  delta = tp->mean[fi] - mean;
	  mean += delta/(double)meann;

	}
      }

      zoffset[fi] = mean;
    }
  }

  bool single_frequency_predictions(size_t frequency_index,
				    const double *model,
				    double offset)
  {
    if (!upscale(model, offset, winverse, super_image)) {
      return false;
    }

    if (linear) {

      for (auto t : traces) {

	for (auto tp : t->paths) {

	  double ttime = tp->linearweights.evaluate_velocityfield(super_image);
	  if (traveltime) {
	    tp->pred[frequency_index] = ttime;
	  } else {
	    tp->pred[frequency_index] = tp->distkm/ttime;
	  }
	  
	}
	
      }
      
    } else {
      
      velocity->set_field(super_image);
      
      for (auto t : traces) {
	
	t->traveltime.construct_traveltime_field();
	
	for (auto tp : t->paths) {
	  
	  //
	  // Get the traveltime
	  //
	  if (traveltime) {
	    tp->pred[frequency_index] = t->traveltime.get_traveltime(tp->coord);
	  } else {
	    tp->pred[frequency_index] = tp->distkm/t->traveltime.get_traveltime(tp->coord);
	  }
	}
      }
    }

    return true;
  }

  double single_frequency_likelihood(size_t frequency_index,
				     const double *model,
				     double offset, 
				     const hierarchicalmodel *hmodel,
				     double *residuals,
				     double *residuals_normed,
				     double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;

    sum = 0.0;

    //
    // First compute predictions
    //
    if (!single_frequency_predictions(frequency_index, model, offset)) {
      return LARGE_LIKELIHOOD;
    }

    //
    // Second loop for computing residuals and likelihood
    //
    for (auto t: traces) {
      for (auto tp: t->paths) {

	if (traveltime) {
	  *res = tp->tt_obs[frequency_index] - tp->pred[frequency_index];
	  
	  
	  sum += hmodel->nll(res,
			     &tp->tt_stddev[frequency_index],
			     1,
			     resn,
			     log_normalization);

	} else {
	  *res = tp->mean[frequency_index] - tp->pred[frequency_index];
	  
	  
	  sum += hmodel->nll(res,
			     &tp->stddev[frequency_index],
			     1,
			     resn,
			     log_normalization);
	}
	
	res ++;
	resn ++;
      }
    }

    return sum;
  }

  bool single_frequency_trace_predictions(size_t frequency_index,
					  size_t trace_index,
					  const double *model,
					  double offset)
  {
    if (!upscale(model, offset, winverse, super_image)) {
      return false;
    }
    
    trace *t = traces[trace_index];

    if (linear) {

      for (auto tp : t->paths) {

	double ttime = tp->linearweights.evaluate_velocityfield(super_image);
	if (traveltime) {
	  tp->pred[frequency_index] = ttime;
	} else {
	  tp->pred[frequency_index] = tp->distkm/ttime;
	}
	
      }
      
    } else {
      
      velocity->set_field(super_image);
      
      t->traveltime.construct_traveltime_field();
      
      for (auto tp : t->paths) {
	
	//
	// Get the traveltime
	//
	if (traveltime) {
	  tp->pred[frequency_index] = t->traveltime.get_traveltime(tp->coord);
	} else {
	  tp->pred[frequency_index] = tp->distkm/t->traveltime.get_traveltime(tp->coord);
	}
	
      }
    }
    
    return true;
  }
  
  double single_frequency_trace_likelihood(size_t frequency_index,
					   size_t trace_index,
					   const double *model,
					   double offset,
					   const hierarchicalmodel *hmodel,
					   double *residuals,
					   double *residuals_normed,
					   double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;

    sum = 0.0;
    
    if (!single_frequency_trace_predictions(frequency_index, trace_index, model, offset)) {
      return LARGE_LIKELIHOOD;
    }
    
    //
    // Second loop for computing residuals and likelihood
    //
    trace *t = traces[trace_index];
    for (auto tp: t->paths) {

      *res = tp->mean[frequency_index] - tp->pred[frequency_index];
      
      sum += hmodel->nll(res,
			 &tp->stddev[frequency_index],
			 1,
			 resn,
			 log_normalization);
      
      res ++;
      resn ++;
    }

    return sum;
  }

  double single_frequency_hierarchical_likelihood(size_t frequency_index,
						  const hierarchicalmodel *hmodel,
						  double *residuals,
						  double *residuals_normed,
						  double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;

    sum = 0.0;

    //
    // Just need loop to sum over residuals
    //
    for (auto t: traces) {
      for (auto tp: t->paths) {

	sum += hmodel->nll(res,
			   &tp->stddev[frequency_index],
			   1,
			   resn,
			   log_normalization);
	
	res ++;
	resn ++;
      }
    }

    return sum;
  }
  
  bool single_trace_predictions(size_t trace_index,
				const double *model,
				const double *offset)
  {
    int rowstride = width * height;
    const double *p = model;


    trace *t = traces[trace_index];
    
    if (linear) {

      for (size_t f = 0; f < frequencies.size(); f ++) {
	
	if (!upscale(p, offset[f], winverse, super_image)) {
	  return false;
	}
	
	for (auto tp : t->paths) {
	  
	  double ttime = tp->linearweights.evaluate_velocityfield(super_image);
	  if (traveltime) {
	    tp->pred[f] = ttime;
	  } else {
	    tp->pred[f] = tp->distkm/ttime;
	  }
	  
	}	
	
	p += rowstride;
      }
      
    } else {
      
      for (size_t f = 0; f < frequencies.size(); f ++) {
	
	if (!upscale(p, offset[f], winverse, super_image)) {
	  return false;
	}
	velocity->set_field(super_image);
	
	t->traveltime.construct_traveltime_field();
	
	for (auto tp : t->paths) {
	  
	  //
	  // Get the traveltime
	  //
	  if (traveltime) {
	    tp->pred[f] = t->traveltime.get_traveltime(tp->coord);
	  } else {
	    tp->pred[f] = tp->distkm/t->traveltime.get_traveltime(tp->coord);
	  }
	  
	}
	
	p += rowstride;
      }
    }
    
    return true;
  }
  
  double single_trace_likelihood(size_t trace_index,
				 const double *model,
				 const double *offset,
				 const hierarchicalmodel *hmodel,
				 double *residuals,
				 double *residuals_normed,
				 double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;
    
    //
    // Compute predictions
    //
    if (!single_trace_predictions(trace_index,
				  model,
				  offset)) {
      return LARGE_LIKELIHOOD;
    }
    
    //
    // Second loop for likelihood
    //
    sum = 0.0;
    trace *t = traces[trace_index];
    for (auto tp: t->paths) {
      
      for (size_t f = 0; f < tp->nfreq; f ++) {
	res[f] = tp->mean[f] - tp->pred[f];
      }
      
      sum += hmodel->nll(res,
			 tp->stddev,
			 tp->nfreq,
			 resn,
			 log_normalization);
      
      res += tp->nfreq;
      resn += tp->nfreq;
    }

    return sum;
  }

  bool predictions(const double *model,
		   const double *offset)
  {
    int rowstride = width * height;
    const double *p = model;

    //
    // First loop needed for computing predictions
    //
    if (linear) {
      for (size_t f = 0; f < frequencies.size(); f ++) {
	
	double o = 0.0;
	if (offset != nullptr) {
	  o = offset[f];
	}
	
	if (!upscale(p, o, winverse, super_image)) {
	  return false;
	}

	for (auto t : traces) {
	  
	  for (auto tp : t->paths) {
	    
	    double ttime = tp->linearweights.evaluate_velocityfield(super_image);
	    if (traveltime) {
	      tp->pred[f] = ttime;
	    } else {
	      tp->pred[f] = tp->distkm/ttime;
	    }
	    
	  }
	}
      }
      
    } else {
      for (size_t f = 0; f < frequencies.size(); f ++) {
	
	double o = 0.0;
	if (offset != nullptr) {
	  o = offset[f];
	}
	
	if (!upscale(p, o, winverse, super_image)) {
	  return LARGE_LIKELIHOOD;
	}
	velocity->set_field(super_image);
	
	for (auto t : traces) {
	  
	  t->traveltime.construct_traveltime_field();
	  
	  for (auto tp : t->paths) {
	    
	    //
	    // Get the traveltime
	    //
	    tp->pred[f] = tp->distkm/t->traveltime.get_traveltime(tp->coord);
	    
	  }
	}
	
	p += rowstride;
      }
    }

    return true;
  }    
  
  double likelihood(const double *model,
		    const double *offset,
		    const hierarchicalmodel *hmodel,
		    double *residuals,
		    double *residuals_normed,
		    double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;

    sum = 0.0;

    if (!predictions(model, offset)) {
      return LARGE_LIKELIHOOD;
    }

    //
    // Second loop for computing residuals and likelihood
    //
    for (auto t: traces) {
      for (auto tp: t->paths) {

	for (size_t f = 0; f < tp->nfreq; f ++) {
	  res[f] = tp->mean[f] - tp->pred[f];
	}
	
	sum += hmodel->nll(res,
			   tp->stddev,
			   tp->nfreq,
			   resn,
			   log_normalization);
	
	res += tp->nfreq;
	resn += tp->nfreq;
      }
    }

    return sum;
  }
  
  double hierarchical_likelihood(const double *model,
				 const hierarchicalmodel *hmodel,
				 double *residuals,
				 double *residuals_normed,
				 double &log_normalization)
  {
    double *res = residuals;
    double *resn = residuals_normed;
    double sum;

    sum = 0.0;

    //
    // Just need loop to sum over residuals
    //
    for (auto t: traces) {
      for (auto tp: t->paths) {

	sum += hmodel->nll(res,
			   tp->stddev,
			   tp->nfreq,
			   resn,
			   log_normalization);
	
	res += tp->nfreq;
	resn += tp->nfreq;
      }
    }

    return sum;
  }

  //
  // Used for generating synthetic observations
  //
  void add_frequency(double f)
  {
    frequencies.push_back(f);
  }
    
  void add_station(const char *code, const coordinate &coord)
  {
    if (find_station(code) == nullptr) {
      stations.push_back(station(code, coord));
    }
  }

  bool add_trace(const char *codeA, const char *codeB, bool minimize)
  {
    coordinate cA;
    if (!find_station_coordinate(codeA, cA)) {
      return false;
    }
    
    coordinate cB;
    if (!find_station_coordinate(codeB, cB)) {
      return false;
    }

    double distkm = LonLat<>::distance_km(cA, cB);

    size_t binsize = (stations.size() - 1)/2 + 1;

    trace *traceA = find_trace(codeA);
    if (traceA == nullptr || !minimize || traceA->paths.size() < binsize ) {
    
      if (traceA == nullptr) {
	traceA = new trace(codeA, cA, velocity, xmin, xmax, ymin, ymax, refinement);
	traces.push_back(traceA);
      }
      
      tracepair *tp = traceA->add_tracepair(codeB, cB, frequencies.size());
      tp->set_distance(distkm);
      
      return true;
    }

    
    trace *traceB = find_trace(codeB);
    if (traceB == nullptr || traceB->paths.size() < traceA->paths.size()) {

      if (traceB == nullptr) {
	traceB = new trace(codeB, cB, velocity, xmin, xmax, ymin, ymax, refinement);
	traces.push_back(traceB);
      }

      tracepair *tp = traceB->add_tracepair(codeA, cA, frequencies.size());
      tp->set_distance(distkm);

      return true;
    } else {

      tracepair *tp = traceA->add_tracepair(codeB, cB, frequencies.size());
      tp->set_distance(distkm);

      return true;
    }

    // Unreachable
    return false;
  }

  void dump_traces()
  {
    for (auto t: traces) {
      printf("%s %d :", t->code.c_str(), (int)t->paths.size());

      for (auto tp: t->paths) {
	printf("%s ", tp->code.c_str());
      }

      printf("\n");
    }
  }

  

  void compute_predictions(const double *model, double offset)
  {
    int rowstride = width * height;
    const double *p = model;
    
    for (size_t f = 0; f < frequencies.size(); f ++) {

      if (!upscale(p, offset, winverse, super_image)) {
	throw WAVETOMO2DEXCEPTION("Negative velocity in prediction model");
      }
      velocity->set_field(super_image);

      for (auto t : traces) {

	t->traveltime.construct_traveltime_field();

	for (auto tp : t->paths) {

	  //
	  // Get the traveltime
	  //
	  tp->pred[f] = tp->distkm/t->traveltime.get_traveltime(tp->coord);
	}
      }

      p += rowstride;
    }
  }

  void compute_predictions_linear(const double *model, double offset)
  {
    int rowstride = width * height;
    const double *p = model;

    for (size_t f = 0; f < frequencies.size(); f ++) {

      if (!upscale(p, offset, winverse, super_image)) {
	throw WAVETOMO2DEXCEPTION("Negative velocity in prediction model");
      }

      for (auto t : traces) {

	for (auto tp : t->paths) {

	  double ttime = 0.0;
	  for (auto tpw : tp->linearweights) {

	    ttime += super_image[tpw.index] * tpw.weight;

	  }

	  if (ttime == 0.0) {
	    throw WAVETOMO2DEXCEPTION("0 Travel time");
	  }
	  
	  tp->pred[f] = tp->distkm/ttime;

	}
      }
      
      p += rowstride;
    }
  }

  void copy_predictions()
  {
    for (size_t f = 0; f < frequencies.size(); f ++) {
      
      for (auto t : traces) {
	
	for (auto tp : t->paths) {
	  
	  tp->mean[f] = tp->pred[f];
	  tp->median[f] = tp->pred[f];
	  tp->mode[f] = tp->pred[f];
	  
	}
      }
    }
  }
  
  void add_noise(Rng &random, hierarchicalmodel &model)
  {
    for (size_t f = 0; f < frequencies.size(); f ++) {

      double sigma = model.noise();
      
      for (auto t : traces) {

	for (auto tp : t->paths) {

	  
	  tp->mean[f] += random.normal(sigma);
	  tp->stddev[f] = sigma;
	  
	}
      }
    }
  }

  bool save_slice_residuals(const char *filename,
			    int n,
			    const double *residuals,
			    const double *normed_residuals)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return -1;
    }

    int i = 0;
    for (auto t : traces) {
      
      for (auto tp : t->paths) {

	if (i >= n) {
	  return false;
	}
	
	fprintf(fp, "%3d %s %s %15.9f %15.9f %15.9f\n",
		i,
		t->code.c_str(),
		tp->code.c_str(),
		tp->distkm,
		residuals[i],
		normed_residuals[i]);
	
	i ++;
      }
    }

    fclose(fp);
    return true;
  }

  bool save_slice_predictions(const char *filename, int slice)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return -1;
    }

    int i = 0;
    for (auto t : traces) {
      
      for (auto tp : t->paths) {

	fprintf(fp, "%3d %s %s %15.9f %15.9f %15.9f %15.9f\n",
		i,
		t->code.c_str(),
		tp->code.c_str(),
		tp->distkm,
		tp->pred[slice],
		tp->mean[slice],
		tp->stddev[slice]);
	
	i ++;
      }
    }

    fclose(fp);
    return true;
  }

  bool upscale(const double *image, double offset, generic_lift_inverse1d_step_t winverse, double *super_image)
  {
    //
    // Copy image to large buffer
    //
    
    memset(super_image, 0, sizeof(double) * super_size);

    double *p = super_image;
    const double *q = image;
    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	p[i] = q[i];
      }

      p += super_width;
      q += width;
    }

    //
    // 2D Inverse wavelet transform steps
    //
    int wwidth = width;
    int wheight = height;

    for (int i = 0; i < super_resolution; i ++) {
      wwidth *= 2;
      wheight *= 2;

      generic_lift_inverse2d_step(super_image,
				  wwidth, wheight,
				  super_width,
				  super_workspace,
				  winverse, winverse);
    }

    for (int i = 0; i < super_size; i ++) {
      super_image[i] += offset;
      
      if (super_image[i] < 0.0) {
	//printf("Negative here\n");
	return false;
      }
    }
    

    return true;
  }

  bool save_super_image(const char *filename) const
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }
    
    for (int j = 0; j < super_height; j ++) {
      for (int i = 0; i < super_width; i ++) {

	fprintf(fp, "%10.6f ", super_image[j * super_width + i]);

      }
      fprintf(fp, "\n");
    }

    fclose(fp);
    return true;
  }

  bool compute_linear_weights()
  {
    for (auto t : traces) {

      // Point A is t->coord
      
      for (auto tp : t->paths) {

	// Point B is tp->coord
	if (!tp->linearweights.compute_weights(xmin, xmax,
					       ymin, ymax,
					       super_width, super_height,
					       RADIUS,
					       t->coord.X(), t->coord.Y(),
					       tp->coord.X(), tp->coord.Y())) {
	  return false;
	}

	// int c = tp->linearweights.ncells();
	// double d = tp->linearweights.total_weight();
	// printf("Linear: %s %s %15.6f %15.6f %2d: %15.9f\n", t->code.c_str(), tp->code.c_str(), tp->distkm, d, c, tp->distkm - d);
      }
    }

    return true;
  }

  bool save_weight_image(const char *filename)
  {
    //
    // Allocate and initialise weight image
    //
    int size = width * height;
    double *weight = new double[size];

    for (int i = 0; i < size; i ++) {
      weight[i] = 0.0;
    }

    //
    // Paint weights into image
    //
    for (auto t : traces) {
      for (auto tp : t->paths)  {

	std::vector<int> indices;
	std::vector<double> weights;
	
	tp->linearweights.get(indices, weights);
	
	for (int i = 0; i < (int)indices.size(); i ++) {
	  weight[indices[i]] += weights[i];
	}
	
      }
    }
      
    //
    // Save image
    //
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "save_weight_image: failed to create file\n");
      return false;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {

	fprintf(fp, "%15.9f ", weight[j * width + i]);
      }

      fprintf(fp, "\n");
    }

    fclose(fp);
    delete weight;

    return true;
  }

 
  bool save_hitcount_image(const char *filename)
  {
    //
    // Allocate and initialise hit count image
    //
    int size = width * height;
    int *hc = new int[size];

    for (int i = 0; i < size; i ++) {
      hc[i] = 0;
    }

    //
    // Paint weights into hit count image
    //
    for (auto t : traces) {
      for (auto tp : t->paths)  {

	std::vector<int> indices;
	std::vector<double> weights;
	
	tp->linearweights.get(indices, weights);
	
	for (int i = 0; i < (int)indices.size(); i ++) {
	  hc[indices[i]] ++;
	}
	
      }
    }
      
    //
    // Save image
    //
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "save_hitcount_image: failed to create file\n");
      return false;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {

	fprintf(fp, "%d ", hc[j * width + i]);
      }

      fprintf(fp, "\n");
    }

    fclose(fp);
    delete hc;

    return true;
  }
  
  bool save_identity_projection_matrix(const char *filename)
  {
    //
    // Allocate and initialise weight image
    //
    int size = width * height;
    double *weight = new double[size];

    for (int i = 0; i < size; i ++) {
      weight[i] = 0.0;
    }

    //
    // Save matrix
    //
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "save_identity_projection_matrix: failed to create file\n");
      return false;
    }
    
    //
    // Paint weights into image
    //
    for (auto t : traces) {
      for (auto tp : t->paths)  {

	std::vector<int> indices;
	std::vector<double> weights;
	
	tp->linearweights.get(indices, weights);

	for (int i = 0; i < size; i ++) {
	  weight[i] = 0.0;
	}
	
	for (int i = 0; i < (int)indices.size(); i ++) {
	  weight[indices[i]] += weights[i];
	}

	for (int i = 0; i < size; i ++) {
	  fprintf(fp, "%16.9e ", weight[i]);
	}
	fprintf(fp, "\n");
      }
    }
      
    fclose(fp);
    delete weight;

    return true;
  }

  bool save_b_dist_error_matrix(const char *filename, int slice)
  {
    //
    // Save vectors together
    //
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "save_identity_projection_matrix: failed to create file\n");
      return false;
    }
    
    //
    // Paint weights into image
    //
    for (auto t : traces) {
      for (auto tp : t->paths)  {

	fprintf(fp, "%16.9e %16.9e %16.9e\n",
		tp->mean[slice], tp->distkm, tp->stddev[slice]);
	
      }
    }
      
    fclose(fp);

    return true;
  }

  bool save_residuals(size_t frequency_index,
		      const char *filename,
		      const double *residuals,
		      const double *residuals_normed)
  {
    const double *res = residuals;
    const double *resn = residuals_normed;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }
    

    for (auto t: traces) {
      for (auto tp: t->paths) {

	double v = tp->mean[frequency_index];
	double tt = tp->distkm/v;
	double deltat = tp->distkm/(v - *res) - tt;

	fprintf(fp, "%s %s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n",
		t->code.c_str(),
		tp->code.c_str(),
		tp->distkm,
		*res,
		v,
		*resn,
		deltat,
		tt);
	
	res ++;
	resn ++;
	
      }
    }

    fclose(fp);
    return true;
  }

  double get_xmin() const {
    return xmin;
  }
  
  double get_xmax() const {
    return xmax;
  }

  double get_ymin() const {
    return ymin;
  }
  
  double get_ymax() const {
    return ymax;
  }

  double get_station_xmin() const {
    return station_xmin;
  }
  
  double get_station_xmax() const {
    return station_xmax;
  }

  double get_station_ymin() const {
    return station_ymin;
  }
  
  double get_station_ymax() const {
    return station_ymax;
  }
  
private:

  class station;
  class trace;
  class tracepair;

  station *find_station(const char *code) 
  {
    for (auto &s : stations) {
      if (s.code == code) {
	return &s;
      }
    }

    return nullptr;
  }
    
  bool find_station_coordinate(const char *code, coordinate &coord)
  {
    station *s = find_station(code);
    if (s != nullptr) {
      coord = s->location;
      return true;
    }

    return false;
  }

  trace *find_trace(const char *code)
  {
    for (auto t: traces) {
      if (t->code == code) {
	return t;
      }
    }

    return nullptr;
  }

  tracepair *find_tracepair(const char *codeA, const char *codeB)
  {
    for (auto t : traces) {

      if (t->code == codeA) {

	tracepair *tp = t->find_tracepair(codeB);
	if (tp != nullptr) {
	  return tp;
	}
	
      }

      if (t->code == codeB) {

	tracepair *tp = t->find_tracepair(codeA);
	if (tp != nullptr) {
	  return tp;
	}

      }

    }

    return nullptr;
  }

  class station {
  public:
    station(const char *_code,
	    const coordinate &_location) :
      code(_code),
      location(_location)
    {
    }
    
    std::string code;
    coordinate location;
  };


  std::vector<station> stations;
  std::vector<double> frequencies;

  class tracepair {
  public:
    tracepair(const char *_code,
	      const coordinate &_coord,
	      double _distkm,
	      size_t _nfreq) :
      code(_code),
      coord(_coord),
      distkm(_distkm),
      nfreq(_nfreq)
    {
      mean = new double[nfreq];
      median = new double[nfreq];
      mode = new double[nfreq];
      stddev = new double[nfreq];
      pred = new double[nfreq];

      tt_stddev = new double[nfreq];
      tt_obs = new double[nfreq];
    }

    ~tracepair()
    {
      delete [] mean;
      delete [] median;
      delete [] mode;
      delete [] stddev;
      delete [] pred;

      delete [] tt_stddev;
      delete [] tt_obs;
    }

    void set_distance(double _distkm)
    {
      distkm = _distkm;
    }

    void set_frequency_observation(int i,
				   double _mean,
				   double _median,
				   double _mode,
				   double _stddev)
    {
      mean[i] = _mean;
      median[i] = _median;
      mode[i] = _mode;
      stddev[i] = _stddev;
    }

    std::string code;
    coordinate coord;
    double distkm;

    size_t nfreq;
    double *mean;
    double *median;
    double *mode;
    double *stddev;

    double *tt_stddev;
    double *tt_obs;
    
    double *pred;

    LinearWeights linearweights;
  };
    
  class trace {
  public:
    trace(const char *_code,
	  const coordinate &_coord,
	  VelocityField<double> *velocity,
	  double xmin,
	  double xmax,
	  double ymin,
	  double ymax,
	  int refinement) :
      code(_code),
      coord(_coord),
      traveltime(coordinate(xmin, ymin),
		 coordinate(xmax, ymax),
		 _coord,
		 *velocity,
		 refinement)
    {
    }

    ~trace()
    {
      for (auto tp: paths) {
	delete tp;
      }
    }

    tracepair* add_tracepair(const char *code,
			     const coordinate &coord,
			     int nfreq)
    {
      tracepair *tp = new tracepair(code, coord, -1.0, nfreq);

      paths.push_back(tp);

      return tp;
    }

    tracepair* find_tracepair(const char *code)
    {
      for (auto tp : paths) {
	if (tp->code == code) {
	  return tp;
	}
      }

      return nullptr;
    }

    std::string code;
    coordinate coord;

    TravelTimeField<LonLat<>, double> traveltime;

    std::vector<tracepair *> paths;
  };

  std::vector<trace*> traces;

  int width;
  int height;
  int super_resolution;

  int super_width;
  int super_height;
  int super_size;
  double *super_image;
  double *super_workspace;

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  
  double station_xmin;
  double station_xmax;
  double station_ymin;
  double station_ymax;

  VelocityField<double> *velocity;
  int refinement;

  bool traveltime;
  bool linear;

  generic_lift_inverse1d_step_t winverse;
  
};

#endif // wavetomo2dobservations_hpp
