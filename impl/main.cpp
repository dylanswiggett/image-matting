#include <iostream>
#include "GridSolver.hpp"
#include "GridSampler.hpp"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

#define SIZE 100

int main(int argc, char** argv) {
  SparseMatrix<double> *m;

  m = new SparseMatrix<double>(SIZE,SIZE);
  m->reserve(SIZE*1000);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (i == j)
        m->insert(j,i) = 2;
      else if (i - j == 1 || i - j == -1)
        m->insert(j,i) = 1;
    }
  }
  m->makeCompressed();

  // std::cout << "'A' Matrix:\n--------------------------------\n" << *m << "\n--------------------------------" << std::endl;

  SparseMatrix<double> *expected;

  expected = new SparseMatrix<double>(SIZE,1);

  for (int i = 0; i < SIZE; i++)
    (*expected).insert(i,0) = 0;

  GridSolver solver(m,expected);

  SparseMatrix<double> *guess = new SparseMatrix<double>(SIZE,1);

  for (int i = 0; i < SIZE; i++)
    (*guess).insert(i,0) = 1000;

  std::cout << "Stepping once... ";

  solver.solve(guess, 1);

  std::cout << "Done!\nexpected has error of:" << std::endl;

  std::cout << (*m) * (*guess) - *expected << std::endl;
}
