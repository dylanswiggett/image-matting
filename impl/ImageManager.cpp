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
  throw "GetLaplacian not implemented";
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleMatrix() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w,image->h);

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      m->insert(x,y) = ((double)GetPixel(x,y)) / 255.0;
    }
  }

  return m;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleVector() {

}

void ImageManager::SaveTo(std::string path) {
  throw "SaveTo not implemented";
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
        return p[0];

    case 4:
        return *(Uint32 *)p;
        break;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}
