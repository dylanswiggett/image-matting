#ifndef _GRID_SAMPLER_HPP_
#define _GRID_SAMPLER_HPP_

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

class GridSampler {
 public:
  GridSampler(const SparseMatrix<double>* mat) : mat_(mat) {};
  SparseMatrix<double>* upsample(int newWidth, int newHeight);
  SparseMatrix<double>* downsample();
 private:
  void upsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int i, int j);
  void downsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int i, int j);
  const SparseMatrix<double>* mat_;
};

#endif  // _GRID_SAMPLER_HPP_
