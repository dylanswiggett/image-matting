#ifndef IMAGE_MANAGER_HPP_
#define IMAGE_MANAGER_HPP_

#include <string>
#include <eigen3/Eigen/Sparse>
#include "SDL/SDL.h"

using namespace Eigen;

class ImageManager {
 public:
  ImageManager(std::string path);
  ~ImageManager();

  SparseMatrix<double,RowMajor>* GetLaplacian();
  SparseMatrix<double,RowMajor>* GetGreyscaleMatrix();
  SparseMatrix<double,RowMajor>* GetGreyscaleVector();

  void SaveTo(std::string path);

 private:
  int GetPixel(int x, int y);

  SDL_Surface *image;
};

#endif  // IMAGE_MANAGER_HPP_