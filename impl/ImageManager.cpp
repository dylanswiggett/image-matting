#include "ImageManager.hpp"

#include "SDL/SDL.h"
#include <iostream>
#include <string>
#include <cmath> 
#include <vector>

#define LAPLACIAN_RAD 1
#define MAX_LAP_RAD (double)(2 * LAPLACIAN_RAD + 1)
#define W_K (MAX_LAP_RAD * MAX_LAP_RAD)
#define EPSILON .01

using namespace Eigen;


ImageManager::ImageManager(std::string path) {
  image = SDL_LoadBMP(path.c_str());
}

ImageManager::ImageManager(SparseMatrix<double,RowMajor>* basis) {
  image = SDL_CreateRGBSurface(0,             // flags
                               basis->rows(), // width
                               basis->cols(), // height
                               32,            // depth
                               0xff000000,    // Rmask
                               0x00ff0000,    // Gmask
                               0x0000ff00,    // Bmask
                               0x000000ff);   // Amask

  Uint32 *pixels = (Uint32 *)image->pixels;

  for (int x = 0; x < basis->rows(); x++) {
    for (int y = 0; y < basis->cols(); y++) {
      int color = (int) (basis->coeff(x,y) * 255.0);
      if (color < 0)
        color = -color;
      if (color > 255)
        color = 255;
      pixels[(y * image->w) + x] = (color << 8) + 0xff;
      // pixels[(y * image->w) + x] = 0xffffffff;
    }
  }
}

ImageManager::~ImageManager() {
  SDL_FreeSurface(image);}

SparseMatrix<double,RowMajor>* ImageManager::GetLaplacian() {
  SparseMatrix<double,RowMajor> *L = new SparseMatrix<double,RowMajor>(image->w * image->h, image->w * image->h);
  // L->reserve(image->w * image->h * 100);

  std::vector<Triplet<double>> tripletList;
  tripletList.reserve(image->w * image->h * 100);

  int centerIndex;

  for (int x1 = 0; x1 < image->w; ++x1) {
    std::cout << x1 << std::endl;
    for (int y1 = 0; y1 < image->h; ++y1) {
      int minx = fmax(x1 - MAX_LAP_RAD, 0);
      int maxx = fmin(x1 + MAX_LAP_RAD + 1, image->w);
      double diagonal = 0;
      for (int x2 = minx; x2 < maxx; ++x2) {
        for (int y2 = 0; y2 < image->h; ++y2) {
          double Lval = LaplaciantAt(x1, y1, x2, y2);
          if (x1 == x2 && y1 == y2)
            centerIndex = tripletList.size();
          else
            diagonal += Lval;
          
          if (Lval != 0 || (x1 == x2 && y1 == y2)) {
            tripletList.push_back(Triplet<double>(x1 * image->h + y1, x2 * image->h + y2, Lval));
          }
        }
      }
      // double diagonal = 0;
      // for (int x2 = 0; x2 < image->w; ++x2) {
      //   for (int y2 = 0; y2 < image->h; ++y2) {
      //     if (x2 != x1 || y2 != y1)
      //       diagonal += L->coeff(x1 * image->h + y1, x2 * image->h + y2);
      //   }
      // }

      tripletList[centerIndex] = Triplet<double>(x1 * image->h + y1, x1 * image->h + y1, -diagonal);
    }
  }

  L->setFromTriplets(tripletList.begin(), tripletList.end());

  // for (int x1 = 0; x1 < image->w; ++x1) {
  //   std::cout << x1 << std::endl;
  //   for (int y1 = 0; y1 < image->h; ++y1) {
  //     double diagonal = 0;
  //     for (int x2 = 0; x2 < image->w; ++x2) {
  //       for (int y2 = 0; y2 < image->h; ++y2) {
  //         if (x2 != x1 || y2 != y1)
  //           diagonal += L->coeff(x1 * image->h + y1, x2 * image->h + y2);
  //       }
  //     }
  //     L->coeffRef(x1 * image->h + y1, x1 * image->h + y1) = -diagonal;
  //   }
  // }

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
  SDL_SaveBMP(image, path.c_str());
}

double ImageManager::LaplaciantAt(int x1, int y1, int x2, int y2) {
  // Based on "A Closed Form Solution to Natural Image Matting" by Levin et al.

  // if (x1 == 0 && y1 == 0)
  //   std::cout << "x1: " << x1 << ", y1: " << y1 << ", x2: " << x2 << ", y2: " << y2 << std::endl;

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

  if (x2 - x1 > MAX_LAP_RAD || y2 - y1 > MAX_LAP_RAD)
    return 0;
  // return 1;

  double I1 = GetIntensity(x1, y1);
  double I2 = GetIntensity(x2, y2);

  double kronecker = (x1 == x2 && y1 == y2) ? 1 : 0;

  double q = 0;

  for (int kx = x2 - LAPLACIAN_RAD; kx <= x1 + LAPLACIAN_RAD; ++kx) {
    for (int ky = y2 - LAPLACIAN_RAD; ky <= y1 + LAPLACIAN_RAD; ++ky) {

      double wk = 0;

      double mean = 0;
      double variance = 0;

      // A single window (centered at pixel k) that contains both pixels i and j
      for (int x = -LAPLACIAN_RAD; x <= LAPLACIAN_RAD; ++x) {
        for (int y = -LAPLACIAN_RAD; y <= LAPLACIAN_RAD; ++y) {
          int xp = kx + x;
          int yp = ky + y;
          if (xp >= 0 && yp >= 0 && xp < image->w && yp < image->h) {
            double intensity = GetIntensity(xp, yp);
            mean += intensity;
            variance += intensity * intensity;
            wk++;
          }
        }
      }

      mean /= wk;
      variance = variance / wk - mean * mean;

      double value = (1.0 / wk) *
                     (1.0 + (1.0 / (EPSILON/wk + variance)) *
                                (I1 - mean) * (I2 - mean));

      q += kronecker - value;
    }
  }

  return q;
}

double ImageManager::GetIntensity(int x, int y) {
  // return ((double) GetPixel(x,y)) / 255.0;
  return ((double) ((GetPixel(x,y) & 0x00ff0000) >> 16)) / 255.0;
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
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN)
            return p[0] << 16 | p[1] << 8 | p[2];
        else
            return p[0] | p[1] << 8 | p[2] << 16;
        break;
        // return p[0];  // Just one color component.
        // TODO: Support multi-color images.

    case 4:
        return *(Uint32 *)p;
        break;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}
