#include <iostream>
#include "GridSolver.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

#define SIZE 1000

int main(int argc, char** argv) {
  compressed_matrix<double> *m;

  m = new compressed_matrix<double>(SIZE,SIZE,SIZE*SIZE);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (i == j)
        (*m)(i,j) = SIZE;
      else
        (*m)(i,j) = 1;
      // (*m)(i,j) = (i * SIZE + j + 1);
    }
  }

  vector<double> *v;

  v = new vector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*v)(i) = i;

  std::cout << "Making gridsolver." << std::endl;

  GridSolver solver(m,v);

  vector<double> *guess = new vector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*guess)(i) = 1;

  std::cout << "Stepping once." << std::endl;

  solver.solve(guess, 1);

  std::cout << *guess << std::endl;
}
