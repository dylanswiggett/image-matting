#include <iostream>
#include "GridSolver.hpp"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

#define SIZE 1024

int main(int argc, char** argv) {
  SparseMatrix<double> *m;

  m = new SparseMatrix<double>(SIZE,SIZE);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (i == j)
        (*m).insert(i,j) = SIZE;
      else
        (*m).insert(i,j) = 1;
      // (*m)(i,j) = (i * SIZE + j + 1);
    }
  }

  SparseVector<double> *v;

  v = new SparseVector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*v).insert(i) = 0;

  std::cout << "Making gridsolver." << std::endl;

  GridSolver solver(m,v);

  SparseVector<double> *guess = new SparseVector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*guess).insert(i) = 1000;

  std::cout << "Stepping once." << std::endl;

  solver.solve(guess, 1);

  std::cout << (*m) * (*guess) << std::endl;
}
