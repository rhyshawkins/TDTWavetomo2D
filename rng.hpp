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
#ifndef rng_h
#define rng_h

#include <memory>

//
// A simple wrapper around gsl random number generator
//
class Rng {
public:

  Rng(int seed);
  ~Rng();

  //
  // Integer random
  //
  int uniform(int n);  // k in (0 .. n - 1) uniform
  int jeffreys(int n); // k in (1 .. n) proportional to 1/k

  int select(int nweights, double *weights);
  void shuffle(int nitems, int *items);
  
  //
  // Floating point random
  //
  double uniform();
  double normal(double sigma);
  double gamma(double a, double b);

  //
  // PDF
  //
  static double pdf_normal(double x, double mean, double sigma);

private:

  class impl;
  std::unique_ptr<impl> pimpl;

};

#endif // rng_h
  
