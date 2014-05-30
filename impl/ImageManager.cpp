#include "ImageManager.hpp"

#include "SDL/SDL.h"
#include <iostream>
#include <string>
#include <cmath> 

#define LAPLACIAN_RAD 1
#define W_K ((double) (2 * LAPLACIAN_RAD + 1) * (2 * LAPLACIAN_RAD + 1))
#define EPSILON .01

using namespace Eigen;


ImageManager::ImageManager(std::string path) {
  image = SDL_LoadBMP(path.c_str());
}

ImageManager::~ImageManager() {
  SDL_FreeSurface(image);}

SparseMatrix<double,RowMajor>* ImageManager::GetLaplacian() {
  SparseMatrix<double,RowMajor> *L = new SparseMatrix<double,RowMajor>(image->w * image->h, image->w * image->h);

  for (int x1 = 0; x1 < image->w; ++x1) {
    for (int y1 = 0; y1 < image->h; ++y1) {
      for (int x2 = 0; x2 < image->h; ++x2) {
        for (int y2 = 0; y2 < image->h; ++y2) {
          double Lval = LaplaciantAt(x1, y1, x2, y2);
          // if (abs(Lval) > .001)
            L->insert(x1 * image->h + y1, x2 * image->h + y2) = Lval;
        }
      }
    }
  }

  return L;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleMatrix() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w,image->h);

  m->reserve(image->w * image->h);

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      m->insert(x,y) = GetIntensity(x,y);
    }
  }

  return m;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleVector() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w * image->h, 1);

  m->reserve(image->w * image->h);

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      m->insert(x * image->h + y, 0) = GetIntensity(x,y);
    }
  }

  return m;
}

void ImageManager::SaveTo(std::string path) {
  throw "SaveTo not implemented";
}

double ImageManager::LaplaciantAt(int x1, int y1, int x2, int y2) {
  // Based on "A Closed Form Solution to Natural Image Matting" by Levin et al.

  if (x2 < x1) {
    int temp = x2;
    x2 = x1;
    x1 = temp;
  }

  if (y2 < y1) {
    int temp = y2;
    y2 = y1;
    y1 = temp;
  }
  if (x2 - x1 > LAPLACIAN_RAD || y2 - y1 > LAPLACIAN_RAD)
    return 0;
  // return 1;

  double I1 = GetIntensity(x1, y1);
  double I2 = GetIntensity(x2, y2);

  double q = 0;

  for (int kx = x2 - LAPLACIAN_RAD; kx <= x1 + LAPLACIAN_RAD; ++kx) {
    for (int ky = y2 - LAPLACIAN_RAD; ky <= y1 + LAPLACIAN_RAD; ++ky) {

      double mean = 0;
      double variance = 0;

      // A single window (centered at pixel k) that contains both pixels i and j
      for (int x = -LAPLACIAN_RAD; x <= LAPLACIAN_RAD; ++x) {
        for (int y = -LAPLACIAN_RAD; y <= LAPLACIAN_RAD; ++y) {
          double intensity = GetIntensity(kx + x, ky + y);
          mean += intensity;
          variance += intensity * intensity;
        }
      }

      mean /= W_K;
      variance = variance / W_K - mean * mean;

      double kronecker = (x1 == x2 && y1 == y2) ? 1 : 0;

      double value = (1.0 / W_K) *
                     (1.0 + 1.0 / (EPSILON/W_K + variance) *
                                (I1 - mean) * (I2 - mean));

      q += kronecker - value;
    }
  }

  return q;
}

double ImageManager::GetIntensity(int x, int y) {
  return ((double) GetPixel(x,y)) / 255.0;
}

int ImageManager::GetPixel(int x, int y) {
  if (x < 0 || y < 0 || x >= image->w || y >= image->h)
    return 0;

  int bpp = image->format->BytesPerPixel;
  Uint8 *p = (Uint8 *)image->pixels + y * image->pitch + x * bpp;

  switch(bpp) {
    case 1:
        return *p;
        break;

    case 2:
        return *(Uint16 *)p;
        break;

    case 3:
        // if(SDL_BYTEORDER == SDL_BIG_ENDIAN)
        //     return p[0] << 16 | p[1] << 8 | p[2];
        // else
        //     return p[0] | p[1] << 8 | p[2] << 16;
        // break;
        return p[0];  // Just one color component.
        // TODO: Support multi-color images.

    case 4:
        return *(Uint32 *)p;
        break;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}
