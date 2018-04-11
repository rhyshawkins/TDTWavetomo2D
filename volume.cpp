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

#include "volume.hpp"

int
volume_save(const char *filename, const double *volume, int width, int height, int depth)
{
  FILE *fp = fopen(filename, "wb");
  int header[3];
  int slicestride;
  const double *p;
  
  if (fp == NULL) {
    return -1;
  }

  header[0] = width;
  header[1] = height;
  header[2] = depth;
  if (fwrite(header, sizeof(int), 3, fp) != 3) {
    return -1;
  }

  slicestride = width * height;
  p = volume;
  for (int k = 0; k < depth; k ++) {

    if (fwrite(p, sizeof(double), slicestride, fp) != (size_t)slicestride) {
      return -1;
    }

    p += slicestride;
  }

  fclose(fp);
  return 0;
}

double *
volume_load(const char *filename, int &width, int &height, int &depth)
{
  FILE *fp = fopen(filename, "rb");
  int header[3];
  int slicestride;
  double *volume;
  double *p;
  
  if (fp == NULL) {
    return nullptr;
  }

  if (fread(header, sizeof(int), 3, fp) != 3) {
    return nullptr;
  }
  width = header[0];
  height = header[1];
  depth = header[2];

  slicestride = width * height;
  volume = new double[slicestride * depth];
  p = volume;
  
  for (int k = 0; k < depth; k ++) {

    if (fread(p, sizeof(double), slicestride, fp) != (size_t)slicestride) {
      return nullptr;
    }

    p += slicestride;
  }

  fclose(fp);
  return volume;
}

int
volume_save_slice_as_text(const char *filename, const double *volume, int width, int height, int depth, int slice)
{
  int slicestride;
  const double *p;
  
  FILE *fp = fopen(filename, "w");

  if (fp == NULL) {
    return -1;
  }

  slicestride = width * height;
  p = volume + slice * slicestride;

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%15.9f ", p[j * width + i]);

    }

    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

int
volume_save_column_as_text(const char *filename,
			   const double *volume,
			   int width, int height, int depth,
			   int x, int y)
{
  int slicestride;
  
  FILE *fp = fopen(filename, "w");

  if (fp == NULL) {
    return -1;
  }

  slicestride = width * height;

  for (int i = 0; i < depth; i ++) {

    fprintf(fp, "%15.9f\n", volume[i * slicestride + y * width + x]);

  }

  fclose(fp);
  return 0;
}
