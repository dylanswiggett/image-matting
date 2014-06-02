#ifndef IMAGE_MANAGER_HPP_
#define IMAGE_MANAGER_HPP_

#include <string>
#include <eigen3/Eigen/Sparse>
#include "SDL/SDL.h"
#include <map>

using namespace Eigen;

typedef struct Coord {
  Coord(double xp, double yp) : x(xp), y(yp) {}
  double x, y;
  bool operator<(const Coord& other) const {
    if (other.x == x)
      return y < other.y;
    else
      return x < other.x;
  }
} Coord;

class ImageManager {
 public:
  ImageManager(std::string path);
  ImageManager(SparseMatrix<double,RowMajor>* basis);
  ImageManager(SparseMatrix<double,RowMajor>* r, SparseMatrix<double,RowMajor>* g, SparseMatrix<double,RowMajor> *b);
  ~ImageManager();

  SparseMatrix<double,RowMajor>* GetLaplacian();
  SparseMatrix<double,RowMajor>* GetGreyscaleMatrix();
  SparseMatrix<double,RowMajor>* GetGreyscaleVector();

  void SaveTo(std::string path);

  int GetW() { return image->w; }
  int GetH() { return image->h; }

  ImageManager* downsize();

 private:
  double LaplaciantAt(int x1, int y1, int x2, int y2);
  void GetIntensity(int x, int y, double* r, double* g, double* b);
  int GetPixel(int x, int y);

  SDL_Surface *image;
};

#endif  // IMAGE_MANAGER_HPP_