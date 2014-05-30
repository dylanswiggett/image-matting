#include "ImageManager.hpp"

#include "SDL/SDL.h"
#include <string>

using namespace Eigen;

ImageManager::ImageManager(std::string path) {
  image = SDL_LoadBMP(path.c_str());
}

ImageManager::~ImageManager() {
  SDL_FreeSurface(image);
}

SparseMatrix<double,RowMajor>* ImageManager::GetLaplacian() {
  SparseMatrix<double,RowMajor> *L = new SparseMatrix<double,RowMajor>(image->w * image->h, image->w * image->h);

  for (int x1 = 0; x1 < image->w; ++x1) {
    for (int y1 = 0; y1 < image->h; ++y1) {
      for (int x2 = 0; x2 < image->h; ++x2) {
        for (int y2 = 0; y2 < image->h; ++y2) {
          double Lval = LaplaciantAt(x1, y1, x2, y2);
          if (Lval > .001)
            L->insert(x1 * image->w + y1, x2 * image->w + y2) = Lval;
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
      m->insert(x,y) = ((double)GetPixel(x,y)) / 255.0;
    }
  }

  return m;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleVector() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w * image->h, 1);

  m->reserve(image->w * image->h);

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      m->insert(x * image->h + y, 0) = ((double)GetPixel(x,y)) / 255.0;
    }
  }

  return m;
}

void ImageManager::SaveTo(std::string path) {
  throw "SaveTo not implemented";
}

double ImageManager::LaplaciantAt(int x1, int y1, int x2, int y2) {
  // Based on "A Closed Form Solution to Natural Image Matting" by Levin et al.
  return 0;
}

int ImageManager::GetPixel(int x, int y) {
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
