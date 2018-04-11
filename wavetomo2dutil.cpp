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

#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "wavetomo2dutil.hpp"

std::string mkfilename(const char *prefix, const char *file)
{
  if (prefix == nullptr) {
    return std::string(file);
  } else {
    return std::string(prefix) + file;
  }
}

std::string mkfilenamerank(const char *prefix, const char *file, int rank)
{
  char buffer[16];

  sprintf(buffer, "-%03d", rank);
  
  if (prefix == nullptr) {
    return std::string(file) + buffer;
  } else {
    return std::string(prefix) + file + buffer;
  }
}

std::string mkformatstring(const char *fmt, ...)
{
  static char *buffer = nullptr;
  static int buffer_size = -1;

  if (buffer == nullptr) {
    buffer_size = 512;
    buffer = new char[buffer_size];
  }

  va_list ap;
  int size;
  
  va_start(ap, fmt);
  size = vsnprintf(buffer, buffer_size, fmt, ap);
  while (size >= buffer_size) {
    delete [] buffer;
    buffer_size *= 2;
    buffer = new char[buffer_size];
    size = vsnprintf(buffer, buffer_size, fmt, ap);
  }
  va_end(ap);

  return std::string(buffer);
}

std::string stripwhitespaceandquotes(const char *s)
{
  int i = 0;
  int j = strlen(s);

  while (i < j && (isspace(s[i]) || s[i] == '"' || s[j] == '\'')) {
    i ++;
  }

  while (j > i && (isspace(s[j]) || s[j] == '"' || s[j] == '\'')) {
    j --;
  }

  return std::string(s + i, j - i - 1);
}

bool loadhierarchicallambda(const char *filename, std::vector<double> &lambda)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return false;
  }

  while (!feof(fp)) {
    double d;
    if (fscanf(fp, "%lf\n", &d) != 1) {
      if (!feof(fp)) {
	return false;
      } else {
	break;
      }
    } else {
      lambda.push_back(d);
    }
  }

  fclose(fp);
  return true;
}

bool save_zoffset(const char *filename,
		  const double *zoffset,
		  int N)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  fprintf(fp, "%d\n", N);

  for (int i = 0; i < N; i ++) {
    fprintf(fp, "%15.9f\n", zoffset[i]);
  }

  fclose(fp);
  return true;
}

double *load_zoffset(const char *filename,
		     int &N)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return nullptr;
  }

  if (fscanf(fp, "%d\n", &N) != 1) {
    return nullptr;
  }

  double *zoffset = new double[N];

  for (int i = 0; i < N; i ++) {

    if (fscanf(fp, "%lf\n", &zoffset[i]) != 1) {
      return nullptr;
    }
  }

  fclose(fp);

  return zoffset;
}

