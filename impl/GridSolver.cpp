#include "./GridSolver.hpp"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>

using namespace boost::numeric::ublas;

GridSolver::GridSolver(const compressed_matrix<double>* aMatrix, const vector<double>* solution) :
  a_mat_(aMatrix), sol_(solution)
{
  // TODO: Check that the matrix and solution have the same dimension. Give appropriate error if not.
}

void GridSolver::solve(vector<double>* guess, double threshold) {
  if (guess->size() != sol_->size())
    throw "Guess wrong size.";
  for (int i = 0; i < 10; i++)
    step(guess);
}

void GridSolver::step(vector<double>* guess) {
  // We use forward substitution to avoid taking inverses.
  // See http://en.wikipedia.org/wiki/Gauss–Seidel_method
  for (size_t i = 0; i < sol_->size(); i++) {
    double sum = (*sol_)(i);
    for (size_t j = 0; j < guess->size(); j++) {
      if (j != i)
        sum -= (*a_mat_)(i,j) * (*guess)(j);
    }
    // sum -= (prod(matrix_row<matrix<double>>(*a_mat_,i),*guess))(0,0);
    (*guess)(i) = (1/(*a_mat_)(i,i)) * sum;
  }
}
