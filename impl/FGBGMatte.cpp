#include "FGBGMatte.hpp"

#include <eigen3/Eigen/Sparse>

#include "ImageManager.hpp"
#include "GridSolver.hpp"

#include <iostream>

using namespace Eigen;

#define GAMMA 10

SparseMatrix<double,RowMajor> *FGBGMatte::GetMatte() {
  // From "Scalable Matting: A Sub-linear Approach," by Lee and Wu.

  SparseMatrix<double,RowMajor> *L = source_->GetLaplacian();
  SparseMatrix<double,RowMajor> *g = guess_->GetGreyscaleVector();

  SparseMatrix<double,RowMajor> f(L->rows(), 1);
  SparseMatrix<double,RowMajor> C(L->rows(), L->cols());

  for (int i = 0; i < g->cols(); i++) {
    double guess = g->coeff(i,0);
    if (guess == 1 || guess == 0) {
      f.insert(i,0) = guess;
      C.insert(i,i) = 1;
    } else {
      f.insert(i,0) = 0;
    }
  }

  f *= GAMMA;

  SparseMatrix<double,RowMajor> A(L->rows(), L->cols());
  A = (*L) + C * GAMMA;

  GridSolver solver(&A, &f);

  for (int i = 0; i < 12; i++) {
    solver.solve(g, 1);
  }

  std::cout << "NORM: " << g->norm() << std::endl;

  // *g *= 1000;

  // for (int i = 0; i < g->rows(); i++) {
  //   if (g->coeff(i,0) < .001) {
  //     g->coeffRef(i,0) = 0;
  //   }
  // }

  return g;
}