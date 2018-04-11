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

#include "rng.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Rng::impl {
public:

  impl(int seed) :
    rng(gsl_rng_alloc(gsl_rng_taus))
  {
    gsl_rng_set(rng, seed);
  }
  
  ~impl()
  {
    gsl_rng_free(rng);
  }

  gsl_rng *rng;
};

Rng::Rng(int seed) :
  pimpl(new impl(seed))
{
}

Rng::~Rng()
{
}

int
Rng::uniform(int n)
{
  return gsl_rng_uniform_int(pimpl->rng, n);
}

int
Rng::jeffreys(int n)
{
  double c;
  c = 0.0;
  for (int i = 1; i <= n; i ++) {
    c += 1.0/(double)i;
  }
  if (c == 0.0) {
    return -1;
  }

  c = 1.0/c;
  double u = uniform();

  int k = 1;
  double pk = c/(double)k;
  while (u > pk) {
    u -= pk;
    k ++;
    pk = c/(double)k;
  }

  return k;
}

int
Rng::select(int nweights, double *weights)
{
  double sum = 0.0;
  for (int i = 0; i < nweights; i ++) {
    sum += weights[i];
  }

  double u = uniform() * sum;

  for (int i = 0; i < nweights; i ++) {
    if (u < weights[i]) {
      return i;
    }

    u -= weights[i];
  }

  return nweights - 1;
}

void
Rng::shuffle(int nitems, int *items)
{
  gsl_ran_shuffle(pimpl->rng, items, nitems, sizeof(int));
}

double
Rng::uniform()
{
  return gsl_rng_uniform(pimpl->rng);
}

double
Rng::normal(double sigma)
{
  return gsl_ran_gaussian_ziggurat(pimpl->rng, sigma);
}

double
Rng::gamma(double a, double b)
{
  return gsl_ran_gamma(pimpl->rng, a, b);
}

double
Rng::pdf_normal(double x, double mean, double sigma)
{
  return gsl_ran_gaussian_pdf(x - mean, sigma);
}
