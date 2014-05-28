#include "./GridSolver.hpp"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

GridSolver::GridSolver(const compressed_matrix<double>* aMatrix, const vector<double>* solution) :
  a_mat_(aMatrix), sol_(solution)
{

}

void GridSolver::solve(const vector<double>* guess, double threshold) {

}
