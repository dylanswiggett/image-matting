#include "ImageManager.hpp"

#include "SDL/SDL.h"
#include <iostream>
#include <string>
#include <cmath> 
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "GridSampler.hpp"

#define LAPLACIAN_RAD 1.0
#define MAX_LAP_RAD (2.0 * LAPLACIAN_RAD + 1.0)
#define W_K (MAX_LAP_RAD * MAX_LAP_RAD)
#define EPSILON .001

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
      pixels[(y * image->w) + x] = (color << 8) + (color << 16) + (color << 24) + 0xff;
    }
  }
}

ImageManager::ImageManager(SparseMatrix<double,RowMajor>* rMat, SparseMatrix<double,RowMajor>* gMat, SparseMatrix<double,RowMajor> *bMat) {
  image = SDL_CreateRGBSurface(0,            // flags
                               rMat->rows(), // width
                               rMat->cols(), // height
                               32,           // depth
                               0,            // Default Rmask
                               0,            // Default Gmask
                               0,            // Default Bmask
                               0);           // Default Amask

  Uint32 *pixels = (Uint32 *)image->pixels;

  for (int x = 0; x < rMat->rows(); x++) {
    for (int y = 0; y < rMat->cols(); y++) {
      int r = (int) (rMat->coeff(x,y) * 255.0);
      int g = (int) (gMat->coeff(x,y) * 255.0);
      int b = (int) (bMat->coeff(x,y) * 255.0);
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;

      Uint32 pixel = SDL_MapRGB(image->format, r, g, b);

      pixels[(y * image->w) + x] = pixel;
    }
  }
}

ImageManager::~ImageManager() {
  SDL_FreeSurface(image);
}

SparseMatrix<double,RowMajor>* ImageManager::GetLaplacian() {
  SparseMatrix<double,RowMajor> *L = new SparseMatrix<double,RowMajor>(image->w * image->h, image->w * image->h);

  std::vector<Triplet<double>> tripletList;
  tripletList.reserve(image->w * image->h * 100);

  for (int x1 = 0; x1 < image->w; ++x1) {
    std::cout << x1 << std::endl;
    for (int y1 = 0; y1 < image->h; ++y1) {
      int minx = fmax(x1 - MAX_LAP_RAD, 0);
      int maxx = fmin(x1 + MAX_LAP_RAD + 1, image->w);
      for (int x2 = minx; x2 < maxx; ++x2) {
        for (int y2 = 0; y2 < image->h; ++y2) {
          double Lval = LaplaciantAt(x1, y1, x2, y2);

          if (Lval != 0) {
            tripletList.push_back(Triplet<double>(x1 * image->h + y1, x2 * image->h + y2, Lval));
          }
        }
      }
    }
  }

  L->setFromTriplets(tripletList.begin(), tripletList.end());

  return L;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleMatrix() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w,image->h);

  std::vector<Triplet<double>> *entries = new std::vector<Triplet<double>>();

  double r, g, b;

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      GetIntensity(x,y,&r,&g,&b);
      entries->push_back(Triplet<double>(x, y, (r + g + b) / 3.0));
    }
  }

  m->setFromTriplets(entries->begin(), entries->end());

  return m;
}

SparseMatrix<double,RowMajor>* ImageManager::GetGreyscaleVector() {
  SparseMatrix<double,RowMajor> *m = new SparseMatrix<double,RowMajor>(image->w * image->h, 1);

  m->reserve(image->w * image->h);

  double r, g, b;

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      GetIntensity(x,y,&r,&g,&b);
      m->insert(x * image->h + y, 0) = (r + g + b) / 3.0;
    }
  }

  return m;
}

void ImageManager::SaveTo(std::string path) {
  SDL_SaveBMP(image, path.c_str());
}

ImageManager* ImageManager::downsize() {
  SparseMatrix<double,RowMajor> *rMat = new SparseMatrix<double,RowMajor>(image->w,image->h);
  SparseMatrix<double,RowMajor> *gMat = new SparseMatrix<double,RowMajor>(image->w,image->h);
  SparseMatrix<double,RowMajor> *bMat = new SparseMatrix<double,RowMajor>(image->w,image->h);

  std::vector<Triplet<double>> *rEntries = new std::vector<Triplet<double>>();
  std::vector<Triplet<double>> *gEntries = new std::vector<Triplet<double>>();
  std::vector<Triplet<double>> *bEntries = new std::vector<Triplet<double>>();

  double r, g, b;

  for (int x = 0; x < image->w; x++) {
    for (int y = 0; y < image->h; y++) {
      GetIntensity(x,y,&r,&g,&b);
      rEntries->push_back(Triplet<double>(x, y, r));
      gEntries->push_back(Triplet<double>(x, y, g));
      bEntries->push_back(Triplet<double>(x, y, b));
    }
  }

  rMat->setFromTriplets(rEntries->begin(), rEntries->end());
  gMat->setFromTriplets(gEntries->begin(), gEntries->end());
  bMat->setFromTriplets(bEntries->begin(), bEntries->end());

  rMat = GridSampler(rMat).downsample();
  gMat = GridSampler(gMat).downsample();
  bMat = GridSampler(bMat).downsample();

  return new ImageManager(rMat, gMat, bMat);
}

double ImageManager::LaplaciantAt(int x1, int y1, int x2, int y2) {
  // Based on "A Closed Form Solution to Natural Image Matting" by Levin et al.

  if (abs(x2 - x1) > MAX_LAP_RAD || abs(y2 - y1) > MAX_LAP_RAD)
    return 0;

  double I1r, I1g, I1b;
  double I2r, I2g, I2b;

  MatrixXd I1(3,1);
  MatrixXd I2(3,1);

  GetIntensity(x1, y1, &I1r, &I1g, &I1b);
  GetIntensity(x2, y2, &I2r, &I2g, &I2b);

  I1 << I1r, I1g, I1b;
  I2 << I2r, I2g, I2b;

  double kronecker = (x1 == x2 && y1 == y2) ? 1 : 0;

  double q = 0;

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

  for (int kx = x2 - LAPLACIAN_RAD; kx <= x1 + LAPLACIAN_RAD; ++kx) {
    for (int ky = y2 - LAPLACIAN_RAD; ky <= y1 + LAPLACIAN_RAD; ++ky) {

      double wk = 0;

      double rm = 0;
      double gm = 0;
      double bm = 0;
      double rrv = 0;
      double rgv = 0;
      double rbv = 0;
      double ggv = 0;
      double gbv = 0;
      double bbv = 0;

      // A single window (centered at pixel k) that contains both pixels i and j
      for (int x = -LAPLACIAN_RAD; x <= LAPLACIAN_RAD; ++x) {
        for (int y = -LAPLACIAN_RAD; y <= LAPLACIAN_RAD; ++y) {
          int xp = kx + x;
          int yp = ky + y;
          if (xp >= 0 && yp >= 0 && xp < image->w && yp < image->h) {
            double r, g, b;
            GetIntensity(xp, yp, &r, &g, &b);

            rm += r;
            gm += g;
            bm += b;

            rrv += r * r;
            rgv += r * g;
            rbv += r * b;
            ggv += g * g;
            gbv += g * b;
            bbv += b * b;
            
            wk++;
          }
        }
      }

      rm  /= wk;
      gm  /= wk;
      bm  /= wk;

      rrv /= wk;
      rgv /= wk;
      rbv /= wk;
      ggv /= wk;
      gbv /= wk;
      bbv /= wk;

      rrv -= rm * rm;
      rgv -= rm * gm;
      rbv -= rm * bm;
      ggv -= gm * gm;
      gbv -= gm * bm;
      bbv -= bm * bm;

      Matrix3d covariance;
      MatrixXd mean(3,1);
      MatrixXd variance(3,1);

      covariance << rrv, rgv, rbv,
                    rgv, ggv, gbv,
                    rbv, gbv, bbv;
      mean << rm, gm, bm;
      variance << rrv, ggv, bbv;

      MatrixXd middleMat = covariance + (EPSILON / wk) * Matrix3d::Identity();

      double value = (1.0 / wk) *
                     (1.0 + ((I1 - mean).transpose() * middleMat.inverse() * (I2 - mean)).coeff(0,0));

      q += kronecker - value;
    }
  }

  return q;
}

void ImageManager::GetIntensity(int x, int y, double *r, double *g, double *b) {
  // return ((double) GetPixel(x,y)) / 255.0;
  // return ((double) ((GetPixel(x,y) & 0x0000ff00) >> 8)) / 255.0;
  Uint8 ri, gi, bi;
  SDL_GetRGB(GetPixel(x,y), image->format, &ri, &bi, &gi);
  *r = ((double) ri) / 255.0;
  *g = ((double) gi) / 255.0;
  *b = ((double) bi) / 255.0;
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

    case 4:
        return *(Uint32 *)p;
        break;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}
