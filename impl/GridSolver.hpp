#ifndef _GRID_SOLVER_HPP_
#define _GRID_SOLVER_HPP_

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

// An interative linear system solver that produces solutions of the form A*u=f
class GridSolver {
 public:
  // Constructs a new grid solver, and performs all heavy calculations necessary
  // for iteration to begin.
  GridSolver(const SparseMatrix<double> *aMatrix, const SparseMatrix<double>* solution);

  // Approximately solves from guess down to some threshold residual.
  // The approximate solution is stored directly back into guess.
  virtual void solve(SparseMatrix<double>* guess, double threshold);
 private:
  void relax(SparseMatrix<double>* guess);

  const SparseMatrix<double> *a_mat_;
  const SparseMatrix<double> *sol_;

  // Gauss-Seidel matrices
  // compressed_matrix<double> p_mat_, g_mat_;
};

#endif  // _GRID_SOLVER_HPP_
