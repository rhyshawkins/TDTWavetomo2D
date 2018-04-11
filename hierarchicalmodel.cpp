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
#include <string.h>

#include "wavetomo2dexception.hpp"

#include "hierarchicalmodel.hpp"

extern "C" {
  #include "slog.h"
};

std::map<std::string, hierarchicalmodel::reader_function_t> hierarchicalmodel::readers =
  {
    {"IndependentGaussian", independentgaussianhierarchicalmodel::read},
  };

hierarchicalmodel::hierarchicalmodel()
{
}

hierarchicalmodel::~hierarchicalmodel()
{
}

hierarchicalmodel *
hierarchicalmodel::load(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("Failed to open file");
    return nullptr;
  }

  char modelname[1024];
  if (fgets(modelname, sizeof(modelname) - 1, fp) == NULL) {
    ERROR("Failed to read hierarchical model name");
    return nullptr;
  }

  // Strip trailing whitespace
  int i = strlen(modelname) - 1;
  while (i > 0 && isspace(modelname[i])) {
    modelname[i] = '\0';
    i --;
  }

  std::map<std::string, hierarchicalmodel::reader_function_t>::iterator r = readers.find(modelname);
  if (r == readers.end()) {
    ERROR("Unknown hierarchical model");
    fclose(fp);
    return nullptr;
  } else {
    hierarchicalmodel *m = (r->second)(fp);
    fclose(fp);

    return m;
  }
}

//
// Independent Gaussian
//
independentgaussianhierarchicalmodel::independentgaussianhierarchicalmodel() :
  lambda(1.0)
{
}

independentgaussianhierarchicalmodel::~independentgaussianhierarchicalmodel()
{
}

int
independentgaussianhierarchicalmodel::nparameters() const
{
  return 1;
}

double
independentgaussianhierarchicalmodel::getparameter(int i) const
{
  if (i == 0) {
    return lambda;
  } else {
    throw WAVETOMO2DEXCEPTION("Invalid index\n");
  }
}

void
independentgaussianhierarchicalmodel::setparameter(int i, double v)
{
  if (i == 0) {
    if (v <= 0.0) {
      throw WAVETOMO2DEXCEPTION("Sigma out of range\n");
    }

    lambda = v;
  } else {
    throw WAVETOMO2DEXCEPTION("Invalid index\n");
  }
}
  
double
independentgaussianhierarchicalmodel::noise() const
{
  return lambda;
}

double
independentgaussianhierarchicalmodel::nll(const double *residuals,
					  const double *sigma,
					  size_t N,
					  double *residuals_normed,
					  double &log_normalization) const
{
  double sum = 0.0;
  double n;
  
  for (size_t i = 0; i < N; i ++) {

    n = sigma[i] * lambda;

    residuals_normed[i] = residuals[i]/n;

    sum += residuals_normed[i] * residuals_normed[i] * 0.5;
    log_normalization += log(n);

  }

  return sum;
}

hierarchicalmodel *
independentgaussianhierarchicalmodel::read(FILE *fp)
{
  double lambda;

  if (fscanf(fp, "%lf", &lambda) != 1) {
    ERROR("Failed to read sigma");
    return nullptr;
  }

  independentgaussianhierarchicalmodel *m = new independentgaussianhierarchicalmodel();
  m->lambda = lambda;

  return m;
}

