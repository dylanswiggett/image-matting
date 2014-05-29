#include "./GridSampler.hpp"

#include <eigen3/Eigen/Sparse>
#include <iostream>

using namespace Eigen;

SparseMatrix<double>* GridSampler::upsample(int newWidth, int newHeight) {
  SparseMatrix<double> *newMat = new SparseMatrix<double>(mat_->rows() * 2, mat_->cols() * 2);

  for (int k = 0; k < mat_->outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(*mat_, k); it; ++it){
      upsampleIndex(mat_, newMat, it.row(), it.col());
    }
  }

  return newMat;
}

SparseMatrix<double>* GridSampler::downsample() {
  int num_rows = mat_->rows() / 2;
  int num_cols = mat_->cols() / 2;
  if (num_rows == 0) num_rows = 1;
  if (num_cols == 0) num_cols = 1;
  SparseMatrix<double> *newMat = new SparseMatrix<double>(num_rows, num_cols);

  for (int k = 0; k < mat_->outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(*mat_, k); it; ++it){
      downsampleIndex(mat_, newMat, it.row(), it.col());
    }
  }

  return newMat;
}

bool sampleValidIndex(int& i, int& j, const SparseMatrix<double>* mat, double *result) {
  if (i >= 0 && i < mat->rows() && j >= 0 && j < mat->cols()) {
    *result += mat->coeff(i,j);
    return true;
  } else {
    return false;
  }
}

void incrValidIndex(SparseMatrix<double>* mat, int i, int j, double amt) {
  if (i >= 0 && i < mat->rows() && j >= 0 && j < mat->cols())
    mat->coeffRef(i,j) += amt;
}

void GridSampler::upsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int iMid, int jMid) {
  int i = iMid * 2;
  int j = jMid * 2;
  int i0 = i - 1;
  int i1 = i + 1;
  int j0 = j - 1;
  int j1 = j + 1;

  double full = src->coeff(iMid, jMid);
  double half = full * .5;
  double quarter = half * .5;

  incrValidIndex(dest, i, j, full);

  incrValidIndex(dest, i0, j0, quarter);
  incrValidIndex(dest, i0, j1, quarter);
  incrValidIndex(dest, i1, j1, quarter);
  incrValidIndex(dest, i1, j0, quarter);

  incrValidIndex(dest, i0, j, half);
  incrValidIndex(dest, i, j0, half);
  incrValidIndex(dest, i1, j, half);
  incrValidIndex(dest, i, j1, half); 
}

void GridSampler::downsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int i, int j) {
  i = (i / 2) * 2;
  j = (j / 2) * 2;
  int i0 = i - 1;
  int i1 = i + 1;
  int j0 = j - 1;
  int j1 = j + 1;

  double result1 = 0;
  double result2 = 0;
  double result3 = 0;
  double scale_factor = 0;

  if (sampleValidIndex(i0,j0, src, &result1)) scale_factor += .0625;
  if (sampleValidIndex(i1,j0, src, &result1)) scale_factor += .0625;
  if (sampleValidIndex(i1,j1, src, &result1)) scale_factor += .0625;
  if (sampleValidIndex(i0,j1, src, &result1)) scale_factor += .0625;

  if (sampleValidIndex(i0,j, src, &result2)) scale_factor += .125;
  if (sampleValidIndex(i1,j, src, &result2)) scale_factor += .125;
  if (sampleValidIndex(i,j0, src, &result2)) scale_factor += .125;
  if (sampleValidIndex(i,j1, src, &result2)) scale_factor += .125;

  if (sampleValidIndex(i,j, src, &result3)) scale_factor += .25;
  dest->coeffRef(i / 2,j / 2) = (result1 * .0625 + result2 * .125 + result3 * .25) / scale_factor;
}
