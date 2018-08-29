#pragma once
#ifndef wavetreemapper_hpp
#define wavetreemapper_hpp

class WavetreeMapper {
public:

  WavetreeMapper(int _kmax)
    kmax(_kmax),
    indices(new int[_kmax]),
    values(new double[_kmax]),
    gradients(new double[_kmax])
  {
  }

  ~WavetreeMapper()
  {
    for (auto &m : mappers) {
      delete m;
    }
  }

  void build_mapper(wavetree_sub_t *tree,
		    int width,
		    int height,
		    generic_lift_inverse1d_step_t xwaveletf,
		    generic_lift_inverse1d_step_t ywaveletf,
		    double *workspace,
		    int index)
  {
    int size = width * height;
    double *image = new double[size];

    memset(image, 0, sizeof(double) * size);

    if (wavetree2d_sub_map_impulse_to_array(tree, index, image, size) < 0) {
      throw WAVETOMO2DEXCEPTION("Failed to create impulse model image\n");
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

    //
    // Check if we can use a bounding box to remove large regions of zeros
    //
    int x0 = 0;
    int x1 = width - 1;
    int y0 = 0;
    int y1 = height - 1;
    bool shrunk;
    
    do {

      shrunk = false;
      
      bool shrinkable = true;
      for (int y = y0; y <= y1; y ++) {
	if (fabs(image[y * width + x0]) > IMAGE_EPSILON) {
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
	if (fabs(image[y * width + x1]) > IMAGE_EPSILON) {
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
	if (fabs(image[y0 * width + x]) > IMAGE_EPSILON) {
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
	if (fabs(image[y1 * width + x]) > IMAGE_EPSILON) {
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
	throw WAVETOMO2DEXCEPTION("Zero impulse response for coeff %d\n", index);
      }
      
    } while (shrunk);

    //
    // Note of caution here, the image allocated above is passed to
    // one of the two classes and is deleted in their destructors
    //
    if (x0 > 0 || x1 < (width - 1) ||
	y0 > 0 || y1 < (height - 1)) {

      return new WavetreeMapBounded(width, height, x0, x1, y0, y1, image);

    } else {

      return new WavetreeMapFull(width, height, image);

    }
  }

  void backproject(wavetree_sub_t *tree,
		   double *dLdI)
  {
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

	p = build_mapper(tree, indices[i]);
	mappers.insert(indices[i], p);

      } else {
	p = j.second;
      }

      gradients[i] = p->backproject(dLdI);
      
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
		    double *_image)
      size(_width * _height),
      image(_image)
    {
    }

    virtual ~WavetreeMapFull()
    {
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

    WavetreeMapFull(int _width,
		    int _height,
		    int _x0, int _x1,
		    int _y0, int _y1,
		    double *_image)
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
      for (y = _y0; y <= _y1; y ++) {
	for (x = _x0; x <= _x1; x ++) {

	  image[(y - _y0) * subwidth + (x - _x0)] = _image[y * width + x];

	}
      }

      delete [] _image;
    }

    virtual ~WavetreeMapFull()
    {
      delete [] image;
    }
    
    virtual double backproject(const double *dLdI) const
    {
      double w = 0.0;
      for (y = y0; y <= y1; y ++) {
	for (x = x0; x <= x1; x ++) {

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
    
  std::map<int, WaveTreeMap *> mappers;

};
