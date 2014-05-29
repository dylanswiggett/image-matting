#include "./GridSampler.hpp"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

SparseMatrix<double>* GridSampler::upsample() {
  SparseMatrix<double> *newMat = new SparseMatrix<double>(mat_->rows() * 2, mat_->cols() * 2);

  for (int k = 0; k < mat_->outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(*mat_, k); it; ++it){
      upsampleIndex(mat_, newMat, it.row(), it.col());
    }
  }

  return newMat;
}

SparseMatrix<double>* GridSampler::downsample() {
  SparseMatrix<double> *newMat = new SparseMatrix<double>(mat_->rows() / 2, mat_->cols() / 2);

  for (int k = 0; k < mat_->outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(*mat_, k); it; ++it){
      downsampleIndex(mat_, newMat, it.row(), it.col());
    }
  }

  return newMat;
}

// void toValidIndex(int& i, int&j, const SparseMatrix<double>* mat) {
//   if (i < 0)
//     i = 0;
//   else if (i >= mat->rows())
//     i = mat->rows() - 1;
//   if (j < 0)
//     j = 0;
//   else if (j >= mat->rows())
//     j = mat->rows() - 1;
// }

void GridSampler::upsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int iMid, int jMid) {
  // for (int iDiff = -1; iDiff <= 1; ++iDiff) {
  //   for (int jDiff = -1; jDiff <= 1; ++jDiff) {
  //     int i = iMid * 2 + iDiff;
  //     int j = jMid * 2 + jDiff;
  //     int i0 = i - 1;
  //     int i1 = i + 1;
  //     int j0 = j - 1;
  //     int j1 = j + 1;
  //     toValidIndex(i0,j0, src);
  //     toValidIndex(i1,j1, src);
  //     toValidIndex(i,j, src);

  //     double result = 0;
  //     result += (src->coeff(i0,j0) + src->coeff(i1,j0) + src->coeff(i1,j1) + src->coeff(i0,j1)) * .25;
  //     result += (src->coeff(i0, j) + src->coeff(i1, j) + src->coeff( i,j0) + src->coeff( i,j1)) * .5;
  //     result += src->coeff(i,j);
  //     dest->coeffRef(i,j) = result;
  //   }
  // }
  int i = iMid * 2;
  int j = jMid * 2;
  int i0 = i - 1;
  int i1 = i + 1;
  int j0 = j - 1;
  int j1 = j + 1;
  toValidIndex(i0,j0, dest);
  toValidIndex(i1,j1, dest);

  double full = src->coeff(iMid, jMid);
  double half = full * .5;
  double quarter = half * .5;

  dest->coeffRef(i, j) += full;
  dest->coeffRef(i0, j0) += quarter; dest->coeffRef(i0, j1) += quarter; dest->coeffRef(i1, j1) += quarter; dest->coeffRef(i1, j0) += quarter;
  dest->coeffRef(i0, j) += half; dest->coeffRef(i, j0) += half; dest->coeffRef(i1, j) += half; dest->coeffRef(i, j1) += half; 
}

void GridSampler::downsampleIndex(const SparseMatrix<double>* src, SparseMatrix<double>* dest, int i, int j) {
  i = (i / 2) * 2;
  j = (j / 2) * 2;
  int i0 = i - 1;
  int i1 = i + 1;
  int j0 = j - 1;
  int j1 = j + 1;
  toValidIndex(i0,j0, src);
  toValidIndex(i1,j1, src);

  double result = 0;
  result += (src->coeff(i0,j0) + src->coeff(i1,j0) + src->coeff(i1,j1) + src->coeff(i0,j1)) * .0625;
  result += (src->coeff(i0, j) + src->coeff(i1, j) + src->coeff( i,j0) + src->coeff( i,j1)) * .125;
  result += src->coeff(i,j) * .25;
  dest->coeffRef(i / 2,j / 2) = result;
}
