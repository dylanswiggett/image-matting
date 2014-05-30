#ifndef _GRID_SOLVER_HPP_
#define _GRID_SOLVER_HPP_

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

// An interative linear system solver that produces solutions of the form A*u=f
class GridSolver {
 public:
  // Constructs a new grid solver, and performs all heavy calculations necessary
  // for iteration to begin.
  GridSolver(const SparseMatrix<double,RowMajor> *aMatrix, const SparseMatrix<double,RowMajor>* solution);

  // Approximately solves from guess down to some threshold residual.
  // The approximate solution is stored directly back into guess.
  virtual void solve(SparseMatrix<double,RowMajor>* guess, double threshold);
 private:
  void relax(SparseMatrix<double,RowMajor>* guess);

  const SparseMatrix<double,RowMajor> *a_mat_;
  const SparseMatrix<double,RowMajor> *sol_;

  // Gauss-Seidel matrices
  // compressed_matrix<double> p_mat_, g_mat_;
};

#endif  // _GRID_SOLVER_HPP_
