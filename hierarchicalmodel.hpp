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
#ifndef hierarchicalmodel_hpp
#define hierarchicalmodel_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <vector>

class hierarchicalmodel {
public:

  hierarchicalmodel();
  virtual ~hierarchicalmodel();

  virtual int nparameters() const = 0;

  virtual double getparameter(int i) const = 0;

  virtual void setparameter(int i, double v) = 0;
  
  virtual double noise() const = 0;
  
  virtual double nll(const double *residual,
		     const double *sigma,
		     size_t N,
		     double *residuals_normed,
		     double &log_normalization) const = 0;

  virtual double nll_gradient(const double *residual,
			      const double *sigma,
			      size_t N,
			      double *residuals_normed,
			      double *weight,
			      double &log_normalization) const = 0;

  static hierarchicalmodel *load(const char *filename);

  typedef hierarchicalmodel* (*reader_function_t)(FILE *fp);

private:
  
  static std::map<std::string, reader_function_t> readers;
  
};

class independentgaussianhierarchicalmodel : public hierarchicalmodel {
public:
  
  independentgaussianhierarchicalmodel();
  virtual ~independentgaussianhierarchicalmodel();

  virtual int nparameters() const;

  virtual double getparameter(int i) const;

  virtual void setparameter(int i, double v);
  
  virtual double noise() const;
    
  virtual double nll(const double *residual,
		     const double *sigma,
		     size_t N,
		     double *residuals_normed,
		     double &log_normalization) const;

  virtual double nll_gradient(const double *residual,
			      const double *sigma,
			      size_t N,
			      double *residuals_normed,
			      double *weight,
			      double &log_normalization) const;

  static hierarchicalmodel *read(FILE *fp);

private:

  double lambda;
  
};

#endif // hierarchicalmodel.hpp
  
