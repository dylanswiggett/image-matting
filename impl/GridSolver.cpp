#include "./GridSolver.hpp"
#include "GridSampler.hpp"

#include <eigen3/Eigen/Sparse>

#include <iostream>

using namespace Eigen;

GridSolver::GridSolver(const SparseMatrix<double> *aMatrix, const SparseMatrix<double>* solution) :
  a_mat_(aMatrix), sol_(solution)
{
  // TODO: Check that the matrix and solution have the same dimension. Give appropriate error if not.
}

void GridSolver::solve(SparseMatrix<double>* guess, double threshold) {
  if (guess->size() != sol_->size())
      throw "Guess wrong size.";


  if (threshold > .01 && guess->rows() > 1) {
    // Relax on Au=f
    step(guess);

    // r <- f - Au
    SparseMatrix<double> residue = *sol_ - *a_mat_ * (*guess);

    SparseMatrix<double> *small_a_mat = GridSampler(a_mat_).downsample();
    SparseMatrix<double> *small_residue = GridSampler(&residue).downsample();

    // e^{2h} = 0
    int downsample_size = small_a_mat->rows();
    SparseMatrix<double> *small_guess = new SparseMatrix<double>(downsample_size, 1);
    for (int i = 0; i < downsample_size; i++)
      (*small_guess).insert(i,0) = 0;

    // e^{2h} = VCYCLE(e^{2h}, r^{2h})
    GridSolver(small_a_mat, small_residue).solve(small_guess, threshold / 2);
    delete small_a_mat;
    delete small_residue;

    // u <- u + e^{h}
    *guess += *GridSampler(small_guess).upsample(guess->rows(), 1);

    // Relax on Au=f
    step(guess);

  }
}

void GridSolver::step(SparseMatrix<double>* guess) {
  // We use forward substitution to avoid taking inverses.
  // See http://en.wikipedia.org/wiki/Gauss–Seidel_method
  for (int i = 0; i < a_mat_->rows(); i++) {
    SparseMatrix<double> result_sum =
      a_mat_->block(0, i, a_mat_->rows(), 1).transpose() * (*guess);
    double sum = sol_->coeff(i,0) - result_sum.coeff(0,0)
               + a_mat_->coeff(i,i) * guess->coeff(i, 0);
    guess->coeffRef(i, 0) = sum / a_mat_->coeff(i,i);
  }
}
