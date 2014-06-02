#ifndef _GRID_SAMPLER_HPP_
#define _GRID_SAMPLER_HPP_

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

class GridSampler {
 public:
  GridSampler(const SparseMatrix<double,RowMajor>* mat) : mat_(mat) {};
  SparseMatrix<double,RowMajor>* upsample(int newWidth, int newHeight);
  SparseMatrix<double,RowMajor>* downsample();
 private:
  double upsampleIndex(const SparseMatrix<double,RowMajor>* src, int i, int j);
  double downsampleIndex(const SparseMatrix<double,RowMajor>* src, int i, int j);
  const SparseMatrix<double,RowMajor>* mat_;
};

#endif  // _GRID_SAMPLER_HPP_
