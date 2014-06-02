#include "./GridSampler.hpp"

#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace Eigen;

SparseMatrix<double,RowMajor>* GridSampler::upsample(int newWidth, int newHeight) {
  SparseMatrix<double,RowMajor> *newMat = new SparseMatrix<double,RowMajor>(newWidth, newHeight);

  std::vector<Triplet<double>> trip_list;
  trip_list.reserve(newWidth * newHeight);

  for (int i = 0; i < newWidth; i++) {
    for (int j = 0; j < newHeight; j++) {
      trip_list.push_back(Triplet<double>(i, j, upsampleIndex(mat_, i, j)));
    }
  }

  newMat->setFromTriplets(trip_list.begin(), trip_list.end());

  return newMat;
}

SparseMatrix<double,RowMajor>* GridSampler::downsample() {
  int num_rows = mat_->rows() / 2;
  int num_cols = mat_->cols() / 2;
  if (mat_->rows() % 2 == 1) num_rows++;
  if (mat_->cols() % 2 == 1) num_cols++;
  SparseMatrix<double,RowMajor> *newMat = new SparseMatrix<double,RowMajor>(num_rows, num_cols);

  std::vector<Triplet<double>> trip_list;
  trip_list.reserve(num_rows * num_cols);

  for (int i = 0; i < newMat->rows(); i++) {
    for (int j = 0; j < newMat->cols(); j++) {
      trip_list.push_back(Triplet<double>(i, j, downsampleIndex(mat_, i, j)));
    }
  }

  newMat->setFromTriplets(trip_list.begin(), trip_list.end());

  return newMat;
}

bool sampleValidIndex(int& i, int& j, const SparseMatrix<double,RowMajor>* mat, double *result) { 
  if (i >= 0 && i < mat->rows() && j >= 0 && j < mat->cols()) {
    *result += mat->coeff(i,j);
    return true;
  } else {
    return false;
  }
}

void incrValidIndex(SparseMatrix<double,RowMajor>* mat, int i, int j, double amt) {
  if (i >= 0 && i < mat->rows() && j >= 0 && j < mat->cols())
    mat->coeffRef(i,j) += amt;
}

double GridSampler::upsampleIndex(const SparseMatrix<double,RowMajor>* src, int iDest, int jDest) {
  int i = (iDest / 2);
  int j = (jDest / 2);
  int i1 = i + 1;
  int j1 = j + 1;

  double result1 = 0;
  double result2 = 0;
  double result3 = 0;
  double scale_factor = 0;

  if (iDest % 2 == 0) {
    if (jDest % 2 == 0) {
      if (sampleValidIndex(i,j, src, &result3)) scale_factor += 1;
    } else {
      if (sampleValidIndex(i,j, src, &result2)) scale_factor += .5;
      if (sampleValidIndex(i,j1, src, &result2)) scale_factor += .5;
    }
  } else {
    if (jDest % 2 == 0) {
      if (sampleValidIndex(i,j, src, &result2)) scale_factor += .5;
      if (sampleValidIndex(i1,j, src, &result2)) scale_factor += .5;
    } else {
      if (sampleValidIndex(i,j, src, &result1)) scale_factor += .25;
      if (sampleValidIndex(i1,j, src, &result1)) scale_factor += .25;
      if (sampleValidIndex(i1,j1, src, &result1)) scale_factor += .25;
      if (sampleValidIndex(i,j1, src, &result1)) scale_factor += .25;
    }
  }

  if (scale_factor == 0)
    return 0;
  else
    return (result1 * .25 + result2 * .5 + result3) / scale_factor;
}

double GridSampler::downsampleIndex(const SparseMatrix<double,RowMajor>* src, int i, int j) {
  i = i * 2;
  j = j * 2;
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

  return (result1 * .0625 + result2 * .125 + result3 * .25) / scale_factor;
}
