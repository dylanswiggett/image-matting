#include <iostream>
#include "GridSolver.hpp"
#include "GridSampler.hpp"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

#define SIZE 10

int main(int argc, char** argv) {
  SparseMatrix<double> *m;

  m = new SparseMatrix<double>(SIZE,SIZE);
  m->reserve(SIZE*1000);
  for (int i = 0; i < SIZE; i++) {
    std::cout << i << std::endl;
    for (int j = 0; j < SIZE; j++) {
      if (i == j)
        m->insert(j,i) = 2;
      else if (i - j == 1 || i - j == -1)
        m->insert(j,i) = 1;
    }
  }
  m->makeCompressed();

  std::cout << *m << std::endl;

  std::cout << *(GridSampler(m).downsample()) << std::endl;

  SparseVector<double> *v;

  v = new SparseVector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*v).insert(i) = i;

  std::cout << "Making gridsolver." << std::endl;

  GridSolver solver(m,v);

  SparseVector<double> *guess = new SparseVector<double>(SIZE);

  for (int i = 0; i < SIZE; i++)
    (*guess).insert(i) = 1000;

  std::cout << "Stepping once." << std::endl;

  solver.solve(guess, 1);

  std::cout << "Done!" << std::endl;

  std::cout << (*m) * (*guess) << std::endl;
}
