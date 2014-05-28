#ifndef _GRID_SOLVER_HPP_
#define _GRID_SOLVER_HPP_

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

// An interative linear system solver that produces solutions of the form A*u=f
class GridSolver {
 public:
  // Constructs a new grid solver, and performs all heavy calculations necessary
  // for iteration to begin.
  GridSolver(const compressed_matrix<double>* aMatrix, const vector<double>* solution);

  // Approximately solves from guess down to some threshold residual.
  // The approximate solution is stored directly back into guess.
  virtual void solve(const vector<double>* guess, double threshold);
 private:
  const compressed_matrix<double> *a_mat_;
  const vector<double> *sol_;

  // Gauss-Seidel matrices
  const compressed_matrix<double> *p_mat_, *g_mat_;
};

#endif  // _GRID_SOLVER_HPP_