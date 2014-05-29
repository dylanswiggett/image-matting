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

  // V-CYCLE
  
  // Relax on Au=f
  // int num_iters = 1000 / guess->rows();
  // if (num_iters == 0) num_iters = 1;
  int num_iters = 1;
  for (int i = 0; i < num_iters; i++)
    relax(guess);

  // std::cout << "INITIAL GUESS WITH THRESH " << threshold << std::endl;
  // std::cout << *guess << std::endl;

  if (threshold > .01 && guess->rows() > 1) {
    // r <- f - Au
    SparseMatrix<double> residue = *sol_ - *a_mat_ * (*guess);

    SparseMatrix<double> *small_a_mat = GridSampler(a_mat_).downsample();
    SparseMatrix<double> *small_residue = GridSampler(&residue).downsample();

    std::cout << "A_MATRIX:\n----------------------\n" << *a_mat_ << "\n----------------------\n" << std::endl;
    std::cout << "RESIDUE:\n----------------------\n" << residue << "\n----------------------\n" << std::endl;
    std::cout << "SMALL A_MATRIX:\n----------------------\n" << *small_a_mat << "\n----------------------\n" << std::endl;
    std::cout << "SMALL RESIDUE:\n----------------------\n" << *small_residue << "\n----------------------\n" << std::endl;

    // e^{2h} = 0
    int downsample_size = small_a_mat->rows();
    SparseMatrix<double> *small_guess = new SparseMatrix<double>(downsample_size, 1);
    for (int i = 0; i < downsample_size; i++)
      (*small_guess).insert(i,0) = 0;

    // e^{2h} = VCYCLE(e^{2h}, r^{2h})
    GridSolver(small_a_mat, small_residue).solve(small_guess, threshold / 2);
    delete small_a_mat;
    delete small_residue;

    std::cout << "FINAL SMALL_ERROR GUESS: " << std::endl;
    std::cout << *small_guess << std::endl;

    // u <- u + e^{h}
    SparseMatrix<double> *error_guess = GridSampler(small_guess).upsample(guess->rows(), 1);
    delete small_guess;

    std::cout << "FINAL ERROR GUESS: " << std::endl;
    std::cout << *error_guess << std::endl;
    *guess += *error_guess;
    delete error_guess;

    std::cout << "FINAL GUESS WITH THRESH " << threshold << std::endl;
    std::cout << *guess << std::endl;
  }

  // Relax on Au=f
  for (int i = 0; i < num_iters; i++)
    relax(guess);
}

void GridSolver::relax(SparseMatrix<double>* guess) {
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
