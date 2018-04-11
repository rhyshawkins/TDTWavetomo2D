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
#ifndef aemutil_hpp

#include <string>
#include <vector>

std::string mkfilename(const char *prefix, const char *file);

std::string mkfilenamerank(const char *prefix, const char *file, int rank);

std::string mkformatstring(const char *fmt, ...);

std::string stripwhitespaceandquotes(const char *s);

bool loadhierarchicallambda(const char *filename, std::vector<double> &lambda);

bool save_zoffset(const char *filename,
		  const double *zoffset,
		  int N);

double *load_zoffset(const char *filename,
		     int &N);

#endif // aemutil.hpp
