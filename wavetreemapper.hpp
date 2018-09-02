#pragma once
#ifndef wavetreemapper_hpp
#define wavetreemapper_hpp

#include <map>

extern "C" {
#include "wavetree2d_sub.h"
};
  
#include "wavetomo2dobservations.hpp"

class WavetreeMapper {
public:

  class WavetreeMap;
  
  static constexpr double IMAGE_EPSILON = 1.0e-9;

  WavetreeMapper(int _kmax) :
    kmax(_kmax),
    indices(new int[_kmax]),
    values(new double[_kmax]),
    gradients(new double[_kmax])
  {
  }

  ~WavetreeMapper()
  {
    for (auto &m : mappers) {
      delete m.second;
    }
  }

  void saveimage(const char *filename, const double *image, int w, int h)
  {
    FILE *fp = fopen(filename, "w");
    for (int j = 0; j < h; j ++) {
      for (int i = 0; i < w; i ++) {
	fprintf(fp, "%16.9e ", image[j * w + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  
  WavetreeMap *build_mapper(wavetree2d_sub_t *tree,
			    int width,
			    int height,
			    generic_lift_inverse1d_step_t xwaveletf,
			    generic_lift_inverse1d_step_t ywaveletf,
			    double *workspace,
			    int index,
			    wavetomo2dobservations<LonLat<>> *obs,
			    double zoffset)
  {
    int size = width * height;
    double *model = new double[size];

    memset(model, 0, sizeof(double) * size);

    if (wavetree2d_sub_map_impulse_to_array(tree, index, model, size) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to create impulse model image\n");
    }

    double maxc = 0.0;
    for (int i = 0; i < size; i ++) {
      if (model[i] > maxc) {
	maxc = model[i];
      }
    }
    

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(model,
			       width,
			       height,
			       width,
			       workspace,
			       xwaveletf,
			       ywaveletf,
			       1 /* SUBTILE always 1 */) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    double maxm = 0.0;
    for (int i = 0; i < size; i ++) {
      if (model[i] > maxm) {
	maxm = model[i];
      }
    }

    //
    // Compute actual image (upscaling)
    //
    int vwidth;
    int vheight;
    const double *vimage = obs->predict_image(model, 0.0, vwidth, vheight);
    if (vimage == nullptr) {
      throw WAVETOMO2DEXCEPTION("Image prediction failed");
    }

    double maxi = 0.0;
    int vsize = vwidth * vheight;
    for (int i = 0; i < vsize; i ++) {
      if (vimage[i] > maxi) {
	maxi = vimage[i];
      }
    }

    //
    // Check if we can use a bounding box to remove large regions of zeros
    //
    int x0 = 0;
    int x1 = vwidth - 1;
    int y0 = 0;
    int y1 = vheight - 1;
    bool shrunk;
    
    do {

      shrunk = false;
      
      bool shrinkable = true;
      for (int y = y0; y <= y1; y ++) {
	if (fabs(vimage[y * vwidth + x0]) > IMAGE_EPSILON) {
	  shrinkable = false;
	  break;
	}
      }
      if (shrinkable) {
	x0 ++;
	shrunk = true;
      }

      shrinkable = true;
      for (int y = y0; y <= y1; y ++) {
	if (fabs(vimage[y * vwidth + x1]) > IMAGE_EPSILON) {
	  shrinkable = false;
	  break;
	}
      }
      if (shrinkable) {
	x1 --;
	shrunk = true;
      }

      shrinkable = true;
      for (int x = x0; x <= x1; x ++) {
	if (fabs(vimage[y0 * vwidth + x]) > IMAGE_EPSILON) {
	  shrinkable = false;
	  break;
	}
      }
      if (shrinkable) {
	y0 ++;
	shrunk = true;
      }
      
      shrinkable = true;
      for (int x = x0; x <= x1; x ++) {
	if (fabs(vimage[y1 * vwidth + x]) > IMAGE_EPSILON) {
	  shrinkable = false;
	  break;
	}
      }
      if (shrinkable) {
	y1 --;
	shrunk = true;
      }

      if (x1 <= x0 ||
	  y1 <= y0) {
	saveimage("model.txt", model, width, height);
	saveimage("upscale.txt", vimage, vwidth, vheight);
	throw WAVETOMO2DEXCEPTION("Zero impulse response for coeff %d %16.9e %16.9e %16.9e\n", index, maxc, maxm, maxi);
      }
      
    } while (shrunk);

    if (x0 > 0 || x1 < (vwidth - 1) ||
	y0 > 0 || y1 < (vheight - 1)) {

      return new WavetreeMapBounded(vwidth, vheight, x0, x1, y0, y1, vimage);

    } else {

      return new WavetreeMapFull(vwidth, vheight, vimage);

    }
  }

  void backproject(wavetree2d_sub_t *tree,
		   int width,
		   int height,
		   generic_lift_inverse1d_step_t xwaveletf,
		   generic_lift_inverse1d_step_t ywaveletf,
		   double *workspace,
		   wavetomo2dobservations<LonLat<>> *obs,
		   double zoffset,
		   double *dLdI)
  {
    int k;
    
    if (wavetree2d_sub_get_model(tree,
				 kmax,
				 indices,
				 values,
				 &k) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to get model values");
    }

    for (int i = 0; i < k; i ++) {

      auto j = mappers.find(indices[i]);
      WavetreeMap *p = nullptr;
      if (j == mappers.end()) {

	p = build_mapper(tree,
			 width, height,
			 xwaveletf, ywaveletf,
			 workspace,
			 indices[i],
			 obs,
			 zoffset);
	mappers.insert(std::pair<int, WavetreeMap*>(indices[i], p));

      } else {
	p = j->second;
      }

      gradients[i] = p->backproject(dLdI);
      printf("Gradient %d: %16.9e\n", i, gradients[i]);
    }
  }

  class WavetreeMap {
  public:

    WavetreeMap()
    {
    }

    virtual ~WavetreeMap()
    {
    }

    virtual double backproject(const double *dLdI) const = 0;
    
  };

  class WavetreeMapFull : public WavetreeMap {
  public:

    WavetreeMapFull(int _width,
		    int _height,
		    const double *_image) :
      size(_width * _height),
      image(new double[_width * _height])
    {
      memcpy(image, _image, size);
    }

    virtual ~WavetreeMapFull()
    {
      delete [] image;
    }
    
    virtual double backproject(const double *dLdI) const
    {
      double w = 0.0;
      for (int i = 0; i < size; i ++) {
	w += dLdI[i] * image[i];
      }
      return w;
    }

    int size;
    double *image;
    
  };

  class WavetreeMapBounded : public WavetreeMap {
  public:

    WavetreeMapBounded(int _width,
		       int _height,
		       int _x0, int _x1,
		       int _y0, int _y1,
		       const double *_image) :
      width(_width),
      height(_height),
      x0(_x0),
      x1(_x1),
      y0(_y0),
      y1(_y1),
      subwidth(_x1 - _x0 + 1),
      subheight(_y1 - _y0 + 1),
      image(nullptr)
    {
      image = new double[subwidth * subheight];
      for (int y = _y0; y <= _y1; y ++) {
	for (int x = _x0; x <= _x1; x ++) {

	  image[(y - _y0) * subwidth + (x - _x0)] = _image[y * width + x];

	}
      }

      delete [] _image;
    }

    virtual ~WavetreeMapBounded()
    {
      delete [] image;
    }
    
    virtual double backproject(const double *dLdI) const
    {
      double w = 0.0;
      for (int y = y0; y <= y1; y ++) {
	for (int x = x0; x <= x1; x ++) {

	  w += image[(y - y0) * subwidth + (x - x0)] * dLdI[y * width + x];

	}
      }

      return w;
    }

    int width, height;
    int x0, x1;
    int y0, y1;
    int subwidth, subheight;
    double *image;
    
  };


  int kmax;
  int *indices;
  double *values;
  double *gradients;
  std::map<int, WavetreeMap *> mappers;

};

#endif // wavetreemapper_hpp
