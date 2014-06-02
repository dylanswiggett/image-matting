#ifndef _FGBGMATTE_HPP_
#define _FGBGMATTE_HPP_

#include "./ImageManager.hpp"
#include <vector>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

class FGBGMatte {
 public:
  FGBGMatte(ImageManager *source, ImageManager *guess);
  ~FGBGMatte();
  SparseMatrix<double,RowMajor> *GetMatte();

 private:
  void buildAMatrix(uint level);
  SparseMatrix<double, RowMajor>* solve(SparseMatrix<double,RowMajor>* guess, SparseMatrix<double,RowMajor>* f, uint level);
  void relax(SparseMatrix<double,RowMajor>* A, SparseMatrix<double,RowMajor>* f, SparseMatrix<double,RowMajor>* guess);

  ImageManager *guess_;
  SparseMatrix<double,RowMajor> *guess_matrix_;

  std::vector<ImageManager *> *image_sizes_;
  std::vector<SparseMatrix<double,RowMajor> *> *a_matrices_;
  std::vector<SparseMatrix<double,RowMajor> *> *f_matrices_;
};

#endif  // _FGBGMATTE_HPP_