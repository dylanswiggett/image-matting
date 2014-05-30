#include <iostream>
#include "GridSolver.hpp"
#include "GridSampler.hpp"
#include "ImageManager.hpp"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

#define SIZE 10000

int main(int argc, char** argv) {
  SparseMatrix<double,RowMajor> *m;

  m = new SparseMatrix<double,RowMajor>(SIZE,SIZE);
  m->reserve(SIZE*3);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (i == j)
        m->insert(i,j) = 2;
      else if (i - j == 1 || i - j == -1)
        m->insert(i,j) = 1;
    }
  }
  m->makeCompressed();

  ImageManager img("test_img/small_2color.bmp");

  std::cout << *(img.GetGreyscaleMatrix()) << std::endl;

  // std::cout << *(GridSampler(GridSampler(m).upsample(7,7)).downsample()) << std::endl;

  return 0;

  // std::cout << "'A' Matrix:\n--------------------------------\n" << *m << "\n--------------------------------" << std::endl;

  SparseMatrix<double,RowMajor> *expected;

  expected = new SparseMatrix<double,RowMajor>(SIZE,1);

  for (int i = 0; i < SIZE; i++)
    (*expected).insert(i,0) = 0;

  GridSolver solver(m,expected);

  SparseMatrix<double,RowMajor> *guess = new SparseMatrix<double,RowMajor>(SIZE,1);

  for (int i = 0; i < SIZE; i++)
    (*guess).insert(i,0) = 1;

  std::cout << "Stepping once... " << std::endl;

  for (int i = 0; i < 4; i++)
    solver.solve(guess, 1);

  std::cout << "Done!\nFound solution:\n" << *guess << "\nWhich has error:" << std::endl;

  std::cout << (*m) * (*guess) - *expected << std::endl;

  delete m;
  delete expected;
  delete guess;
}
