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
#ifndef linearweights_hpp
#define linearweights_hpp

#include <vector>

#include <math.h>

class LinearWeights {
public:

  static constexpr double EPSILON = 1.0e-9;

  double clamp(double x, double target)
  {
    if (fabs(x - target) < EPSILON) {
      return target;
    }
    return x;
  }

  double dclamp(double x, double target1, double target2)
  {
    if (fabs(x - target1) < EPSILON) {
      return target1;
    } else if (fabs(x - target2) < EPSILON) {
      return target2;
    }

    return x;
  }
  
  LinearWeights()
  {
  }

  ~LinearWeights()
  {
  }

  bool compute_weights(double xmin, double xmax,
		       double ymin, double ymax,
		       int width, int height,
		       double radius,
		       double x0, double y0,
		       double x1, double y1)
  {
    //
    // Compute a reasonable threshold distance between points
    //
    double threshold = 1.5 * (xmax - xmin)/(double)width * M_PI/180.0 * radius;

    if (gc_dist(x0, y0, x1, y1, radius) < threshold) {
      
      return add_linesegment(x0, y0,
			     x1, y1,
			     xmin, xmax,
			     ymin, ymax,
			     width, height,
			     radius);

    } else {

      //
      // Split the line in half recursively. This ensures we approximately
      // get a great circle path
      //

      double xc, yc;

      gc_midpoint(x0, y0, x1, y1, xc, yc);
      
      if (!compute_weights(xmin, xmax,
			   ymin, ymax,
			   width, height,
			   radius,
			   x0, y0,
			   xc, yc)) {
	return false;
      }

      return compute_weights(xmin, xmax,
			     ymin, ymax,
			     width, height,
			     radius,
			     xc, yc,
			     x1, y1);
    }    
  }

  bool compute_weights_direct(double xmin, double xmax,
			      double ymin, double ymax,
			      int width, int height,
			      double radius,
			      double x0, double y0,
			      double x1, double y1)
  {
    return add_linesegment(x0, y0,
			   x1, y1,
			   xmin, xmax,
			   ymin, ymax,
			   width, height,
			   radius);
  }

  double evaluate_velocityfield(const double *image) const
  {
    double ttime = 0.0;

    // Weight is distance so travel time equals sum weight[i]/image[index[i]]
    for (auto &w: weights) {

      ttime += w.weight/image[w.idx];
      
    }

    return ttime;
  }

  void clear()
  {
    weights.clear();
  }

  void get(std::vector<int> &_indices, std::vector<double> &_weights)
  {
    for (auto &w: weights) {
      _indices.push_back(w.idx);
      _weights.push_back(w.weight);
    }
  }

  int ncells() const
  {
    return (int)weights.size();
  }
  
  double total_weight() const
  {
    double t = 0.0;
    for (auto &w: weights) {
      t += w.weight;
    }

    return t;
  }

private:

  bool add_linesegment(double x0, double y0,
		       double x1, double y1,
		       double xmin, double xmax,
		       double ymin, double ymax,
		       int width, int height,
		       double radius)
  {
    int c1;
    int c2;
    
    double llon, rlon;
    double ulat, blat;
    
    double lon;
    double lat;
    
    double t;
    
    if (fabs(x1 - x0) < EPSILON &&
	fabs(y1 - y0) < EPSILON) {
      return true;
    }
    
    c1 = cell_index(x0, y0, 
		    xmin,
		    xmax,
		    ymin,
		    ymax,
		    width,
		    height);
    if (c1 < 0) {
      fprintf(stderr, "add_linesegment: point 2 out of range %f %f\n", x0, y0);
      return false;
    }
    
    c2 = cell_index(x1, y1, 
		    xmin,
		    xmax,
		    ymin,
		    ymax,
		    width,
		    height);
    if (c2 < 0) {
      fprintf(stderr, "add_linesegment: point 2 out of range %f %f\n", x1, y1);
      return false;
    }

    if (c1 == c2) {
      if (!add_cell(c1, gc_dist(x0, y0, x1, y1, radius))) {
	return false;
      }

      return true;
    }

    llon = xmin + (double)(c1 % width)*(xmax - xmin)/(double)width;
    rlon = xmin + (double)(1 + c1 % width)*(xmax - xmin)/(double)width;

    if ((x0 == rlon) || (x0 == llon)) {
      /*
       * Deal with an edge case
       */
      if ((c2 % width) > (c1 % width)) {
	return add_linesegment(x0 + EPSILON/2.0, y0,
			       x1, y1, 
			       xmin,
			       xmax,
			       ymin,
			       ymax,
			       width,
			       height,
			       radius);
      } else if ((c2 % width) < (c1 % width)) {
	return add_linesegment(x0 - EPSILON/2.0, y0,
			       x1, y1, 
			       xmin,
			       xmax,
			       ymin,
			       ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    ulat = ymax - (double)((int)(c1/width) )*(ymax - ymin)/(double)height;
    blat = ymax - (double)(1 + (int)(c1/width))*(ymax - ymin)/(double)height;
    
    if ((y0 == ulat) || (y0 == blat)) {
      /*
       * Deal with an edge case
       */
      if ((c2 / width) > (c1 / width)) {
	return add_linesegment(x0, y0 - EPSILON/2.0,
			       x1, y1, 
			       xmin,
			       xmax,
			       ymin,
			       ymax,
			       width,
			       height,
			       radius);
      } else if ((c2 / width) < (c1 / width)) {
	return add_linesegment(x0, y0 + EPSILON/2.0,
			       x1, y1, 
			       xmin,
			       xmax,
			       ymin,
			       ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    /*
     * Check left
     */
    if (x1 < x0) {

      if (fabs(x1 - llon) < EPSILON) {
	t = 1.0;
      } else {
	t = dclamp((llon - x0)/(x1 - x0), 0.0, 1.0);
      }

      if (t > 0.0 && t <= 1.0) {
	
	lat = dclamp(y0 + t * (y1 - y0), ulat, blat);
	
	if (lat <= ulat && lat >= blat) {
	  /*
	   * Exiting cell to left, add the part of the line segment in this cell.
	   */
	  if (!add_cell(c1, gc_dist(x0, y0, llon, lat, radius))) {
	    return false;
	  }
	  
	  if (t < 1.0) {
	    /* Not just touching, need to add rest of the line segment */
	    
	    if (lat == ulat) {
	      lat += EPSILON/2.0;
	    } else if (lat == blat) {
	      lat -= EPSILON/2.0;
	    }
	    
	    if (!add_linesegment(llon - EPSILON/2.0, lat,
				 x1, y1,
				 xmin, xmax,
				 ymin, ymax,
				 width,
				 height,
				 radius)) {
	      return false;
	    }
	  }
	  
	  return true;
	}
      } else if (t == 0.0) {
	return add_linesegment(x0 - EPSILON/2.0, y0,
			       x1, y1,
			       xmin, xmax,
			       ymin, ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    /*
     * Check right
     */
    if (x1 > x0) {

      if (fabs(x1 - rlon) < EPSILON) {
	t = 1.0;
      } else {
	t = dclamp((rlon - x0)/(x1 - x0), 0.0, 1.0);
      }

      if (t > 0.0 && t <= 1.0) {
	
	lat = dclamp(y0 + t * (y1 - y0), ulat, blat);

	if (lat <= ulat && lat >= blat) {

	  /*
	   * Exiting cell to left, add the part of the line segment in this cell.
	   */
	  if (!add_cell(c1, gc_dist(x0, y0, rlon, lat, radius))) {
	    return false;
	  }

	  
	  if (t < 1.0) {
	    /* Not just touching, need to add rest of the line segment */
	    
	    if (lat == ulat) {
	      lat += EPSILON/2.0;
	    } else if (lat == blat) {
	      lat -= EPSILON/2.0;
	    }
            // A
	    // printf("A %15.9f %15.9f %15.9f %15.9f\n", rlon + EPSILON/2.0, lat,
	    // 	   x1, y1);
	    if (!add_linesegment(rlon + EPSILON/2.0, lat,
				 x1, y1,
				 xmin, xmax,
				 ymin, ymax,
				 width,
				 height,
				 radius) ) {
	      fprintf(stderr, "add_linesegment: failed to add right tail\n");
	      return false;
	    }
	  }
	  
	  return true;
	}
      } else if (t == 0.0) {
	return add_linesegment(x0 + EPSILON/2.0, y0,
			       x1, y1,
			       xmin, xmax,
			       ymin, ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    /*
     * Check up
     */
    if (y1 > y0) {

      if (fabs(ulat - y1) < EPSILON) {
	t = 1.0;
      } else {
	t = dclamp((ulat - y0)/(y1 - y0), 0.0, 1.0);
      }

      if (t > 0 && t <= 1.0) {
	
	lon = dclamp(x0 + t * (x1 - x0), rlon, llon);
	
	if (lon <= rlon && lon >= llon) {
	  
	  /*
	   * Exiting cell above, add the part of the line segment in this cell.
	   */
	  if (!add_cell(c1, gc_dist(x0, y0, lon, ulat, radius))) {
	    return false;
	  }
	  
	  if (t < 1.0) {
	    /* Not just touching, need to add rest of the line segment */
	    
	    if (lon == rlon) {
	      lon += EPSILON/2.0;
	    } else if (lon == llon) {
	      lon -= EPSILON/2.0;
	    }
	    // B
	    // printf("B %15.9f %15.9f %15.9f %15.9f\n", lon , ulat + EPSILON/2.0,
	    // 	   x1, y1);
	    if (!add_linesegment(lon, ulat + EPSILON/2.0,
				 x1, y1,
				 xmin, xmax,
				 ymin, ymax,
				 width,
				 height,
				 radius)) {
	      return false;
	    }
	  }
	  
	  return true;
	}
      } else if (t == 0.0) {
	return add_linesegment(x0, y0 + EPSILON/2.0,
			       x1, y1,
			       xmin, xmax,
			       ymin, ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    /*
     * Check down 
     */
    if (y0 > y1) {

      if (fabs(blat - y1) < EPSILON) {
	t = 1.0;
      } else {
	t = dclamp((blat - y0)/(y1 - y0), 0.0, 1.0);
      }
      
      if (t > 0 && t <= 1.0) {
	
	lon = dclamp(x0 + t * (x1 - x0), rlon, llon);
	
	if (lon <= rlon && lon >= llon) {
	  
	  /*
	   * Exiting cell above, add the part of the line segment in this cell.
	   */
	  if (!add_cell(c1, gc_dist(x0, y0, lon, blat, radius))) {
	    return false;
	  }
	  
	  if (t < 1.0) {
	    /* Not just touching, need to add rest of the line segment */

	    if (lon == rlon) {
	      lon += EPSILON/2.0;
	    } else if (lon == llon) {
	      lon -= EPSILON/2.0;
	    }
	    
	    if (!add_linesegment(lon, blat - EPSILON/2.0,
				x1, y1,
				xmin, xmax,
				ymin, ymax,
				width,
				height,
				radius)) {
	      return false;
	    }
	  }
	  
	  return true;
	}
      } else if (t == 0.0) {
	return add_linesegment(x0, y0 - EPSILON/2.0,
			       x1, y1,
			       xmin, xmax,
			       ymin, ymax,
			       width,
			       height,
			       radius);
      }
    }
    
    fprintf(stderr, "add_linesegment: failed to add line segment %.10f %.10f - %.10f %.10f (%d %d)\n",
	    x0, y0,
	    x1, y1,
	    c1, c2);
    fprintf(stderr, "%15.9f %15.9f : %15.9f %15.9f\n", llon, rlon, ulat, blat);
    return false;
  }
  
  int cell_index(double x,
		 double y,
		 double xmin,
		 double xmax,
		 double ymin,
		 double ymax,
		 int width,
		 int height)
  {
    int i;
    int j;
    
    i = (int)floor((x - xmin)/(xmax - xmin) * (double)width);
    j = (int)floor((ymax - y)/(ymax - ymin) * (double)height);
    
    if (i < 0 || i >= width) {
      return -1;
    }
    
    if (j < 0 || j >= height) {
      return -1;
    }
    
    return i + (j * width);
  }

  bool add_cell(int idx,
		double dist)
  {
    for (auto &w: weights) {
      if (w.idx == idx) {
	w.weight += dist;
	return true;
      }
    }

    //
    // Fall through to add a new index/weight pair
    //

    weights.push_back(weight_t(idx, dist));
    
    return true;
  }

  double gc_dist(double x0, 
		 double y0, 
		 double x1, 
		 double y1, 
		 double r)
  {
    double dlon;
    double dlat;
    
    double hsinlat;
    double hsinlon;
    
    x0 *= M_PI/180.0;
    y0 *= M_PI/180.0;
    
    x1 *= M_PI/180.0;
    y1 *= M_PI/180.0;
    
    dlon = x1 - x0;
    dlat = y1 - y0;
    
    hsinlat = sin(dlat/2.0);
    hsinlon = sin(dlon/2.0);
    
    return 2.0 * asin(sqrt(hsinlat * hsinlat + cos(y0)*cos(y1)*(hsinlon*hsinlon))) * r;
  }

  void gc_midpoint(double x0, double y0, double x1, double y1, double &xc, double &yc)
  {
    double phi0 = (90.0 - y0) * M_PI/180.0;
    double theta0 = x0 * M_PI/180.0;

    double phi1 = (90.0 - y1) * M_PI/180.0;
    double theta1 = x1 * M_PI/180.0;
    
    double dtheta = theta1 - theta0;
    double bx = sin(phi1) * cos(dtheta);
    double by = sin(phi1) * sin(dtheta);
    double delta = sin(phi0) + bx;
    
    yc = 180.0/M_PI * atan2(cos(phi0) + cos(phi1), sqrt(delta*delta + by*by));
    xc = 180.0/M_PI * (theta0 + atan2(by, sin(phi0) + bx));
  }
  
  struct weight_t {

    weight_t(int _idx, double _weight) :
      idx(_idx),
      weight(_weight)
    {
    }
      
    int idx;
    double weight;
    
  };
  
  std::vector<weight_t> weights;
  
};

#endif // linearweights_hpp
