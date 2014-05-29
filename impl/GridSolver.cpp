#include "./GridSolver.hpp"

#include <eigen3/Eigen/Sparse>

#include <iostream>

using namespace Eigen;

GridSolver::GridSolver(const SparseMatrix<double> *aMatrix, const SparseVector<double>* solution) :
  a_mat_(aMatrix), sol_(solution)
{
  // TODO: Check that the matrix and solution have the same dimension. Give appropriate error if not.
}

void GridSolver::solve(SparseVector<double>* guess, double threshold) {
  if (guess->size() != sol_->size())
    throw "Guess wrong size.";
  for (int i = 0; i < 1; i++)
    step(guess);
}

void GridSolver::step(SparseVector<double>* guess) {
  // We use forward substitution to avoid taking inverses.
  // See http://en.wikipedia.org/wiki/Gauss–Seidel_method
  for (int i = 0; i < a_mat_->rows(); i++) {
    SparseMatrix<double> result_sum =
      a_mat_->block(0, i, a_mat_->rows(), 1).transpose() * (*guess);
    double sum = sol_->coeff(i) - result_sum.coeff(0,0)
               + a_mat_->coeff(i,i) * guess->coeff(i);
    guess->coeffRef(i) = sum / a_mat_->coeff(i,i);
  }
}
