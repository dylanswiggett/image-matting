#include "FGBGMatte.hpp"

#include <eigen3/Eigen/Sparse>

#include "ImageManager.hpp"
#include "GridSampler.hpp"

#include <iostream>
#include <vector>

using namespace Eigen;

#define GAMMA 10000

FGBGMatte::FGBGMatte(ImageManager *source, ImageManager *guess)
  : guess_(guess)
{
  guess_matrix_ = guess_->GetGreyscaleMatrix();

  image_sizes_ = new std::vector<ImageManager *>();
  a_matrices_ = new std::vector<SparseMatrix<double,RowMajor> *>();
  f_matrices_ = new std::vector<SparseMatrix<double,RowMajor> *>();
  image_sizes_->push_back(source);

  buildAMatrix(0);
}

FGBGMatte::~FGBGMatte() {
  delete a_matrices_;
  delete f_matrices_;
}

SparseMatrix<double,RowMajor> *FGBGMatte::GetMatte() {
  // From "Scalable Matting: A Sub-linear Approach," by Lee and Wu.

  SparseMatrix<double, RowMajor> *guess = guess_->GetGreyscaleMatrix();

  for (int i = 0; i < 30; i++) {
    std::cout << "Iter: " << i << std::endl; //";    Resid. Norm: " << (A * *g).norm() << std::endl;
    guess = solve(guess, (*f_matrices_)[0], 0);
  }

  return guess;
}

void FGBGMatte::buildAMatrix(uint level) {

  std::cout << "Get laplacian" << std::endl;
  SparseMatrix<double,RowMajor> *L = (*image_sizes_)[level]->GetLaplacian();

  SparseMatrix<double,RowMajor> *A = new SparseMatrix<double,RowMajor>();

  SparseMatrix<double,RowMajor> *f = new SparseMatrix<double,RowMajor>((*image_sizes_)[level]->GetW(), (*image_sizes_)[level]->GetH());
  SparseMatrix<double,RowMajor> C(L->rows(), L->cols());

  f->reserve(L->rows());
  C.reserve(L->rows());

  int dim = pow(2,level);

  std::cout << "Build" << std::endl;

  for (int i = 0; i < f->rows(); i++) {
    for (int j = 0; j < f->cols(); j++) {
      double guess_val = guess_matrix_->coeff(i * dim,j * dim);
      if (guess_val >= 1 || guess_val == 0) {
        f->insert(i,j) = guess_val;
        C.insert(i * f->cols() + j, i * f->cols() + j) = 1;
      } else {
        f->insert(i,j) = 0;
      }
    }
  }

  if (level == 0) {
    (*f) *= GAMMA;
    f_matrices_->push_back(f);
  } else {
    delete f;
  }

  (*A) = (*L) + C * GAMMA;

  a_matrices_->push_back(A);

  delete L;

  std::cout << "done" << std::endl;
}

SparseMatrix<double,RowMajor>* matToVec(SparseMatrix<double,RowMajor>* mat, int length) {
  SparseMatrix<double,RowMajor> *vec = new SparseMatrix<double,RowMajor>(length, 1);

  std::vector<Triplet<double>> trip_list;
  trip_list.reserve(length);

  for (int k = 0; k < mat->outerSize(); ++k){
    for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it(*mat, k); it; ++it){
      trip_list.push_back(Triplet<double>(it.row() * mat->cols() + it.col(), 0, it.value()));
    }
  }

  vec->setFromTriplets(trip_list.begin(), trip_list.end());

  return vec;
}

SparseMatrix<double,RowMajor>* vecToMat(SparseMatrix<double,RowMajor>* vec, int rows, int cols) {
  SparseMatrix<double,RowMajor> *mat = new SparseMatrix<double,RowMajor>(rows, cols);

  std::vector<Triplet<double>> trip_list;
  trip_list.reserve(vec->size());

  for (int i = 0; i < vec->rows(); i++){
    int row = i / cols;
    int col = i % cols;
    if (row < rows && col < cols) {
      trip_list.push_back(Triplet<double>(row, col, vec->coeff(i,0)));
    }
  }

  mat->setFromTriplets(trip_list.begin(), trip_list.end());

  return mat;
}

// Must be called with level = 0.
SparseMatrix<double, RowMajor>* FGBGMatte::solve(SparseMatrix<double,RowMajor>* guess, SparseMatrix<double,RowMajor>* f, uint level) {
  SparseMatrix<double,RowMajor>* a_mat;

  // Terminating condition.
  if (f->rows() == 1 || f->cols() == 1)
    return guess;

  // Get our scaled matrices, if necessary.
  if (level >= a_matrices_->size()) {
    // SparseMatrix<double,RowMajor> *downsampled_image = GridSampler((*image_sizes_)[level - 1]->GetGreyscaleMatrix()).downsample();

    image_sizes_->push_back((*image_sizes_)[level - 1]->downsize());
    buildAMatrix(level);
  }
  a_mat = (*a_matrices_)[level];

  SparseMatrix<double,RowMajor> *guess_vec = matToVec(guess, a_mat->rows());
  SparseMatrix<double,RowMajor> *f_vec = matToVec(f, a_mat->rows());

  relax(a_mat, f_vec, guess_vec);

  // r <- f - Au
  SparseMatrix<double,RowMajor> residue = (*f_vec) - (*a_mat) * (*guess_vec);
  std::cout << "Residue: " << residue.norm() << std::endl;
  SparseMatrix<double,RowMajor> *residue_square = vecToMat(&residue, f->rows(), f->cols());
  SparseMatrix<double,RowMajor> *small_residue = GridSampler(residue_square).downsample();

  // e^{2h} = 0
  SparseMatrix<double,RowMajor> *small_guess = new SparseMatrix<double,RowMajor>(small_residue->rows(), small_residue->cols());

  // e^{2h} = VCYCLE(e^{2h}, r^{2h})
  small_guess = solve(small_guess, small_residue, level + 1);

  *small_guess *= 2.0; // Successive over-relaxation.

  // u <- u + e^{h}
  SparseMatrix<double,RowMajor> *error_guess = GridSampler(small_guess).upsample(guess->rows(), guess->cols());
  delete small_guess;

  SparseMatrix<double,RowMajor> final_guess = *guess + *error_guess;
  SparseMatrix<double,RowMajor> *final_guess_vec = matToVec(&final_guess, a_mat->rows());
  delete error_guess;

  // Relax on Au=f
  for (uint i = 0; i <= level; i++)
    relax(a_mat, f_vec, final_guess_vec);

  // Return the guess in image shape again!
  SparseMatrix<double,RowMajor> *final_guess_square = vecToMat(final_guess_vec, guess->rows(), guess->cols());
  return final_guess_square;
}

// Gauss-Seidel Relaxation
void FGBGMatte::relax(SparseMatrix<double,RowMajor>* A, SparseMatrix<double,RowMajor>* f, SparseMatrix<double,RowMajor>* guess) {
  for (int i = 0; i < A->rows(); i++) {
    SparseMatrix<double,RowMajor> row = A->middleRows(i, 1);
    SparseMatrix<double,RowMajor> result = row * (*guess);
    double a_mat_coeff = A->coeff(i,i);
    double sum = f->coeff(i,0) - result.coeff(0,0)
               + a_mat_coeff * guess->coeff(i,0);
    guess->coeffRef(i, 0) = sum / a_mat_coeff;
  }
}
