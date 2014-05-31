#include "FGBGMatte.hpp"

#include <eigen3/Eigen/Sparse>

#include "ImageManager.hpp"
#include "GridSolver.hpp"

#include <iostream>

using namespace Eigen;

#define GAMMA 10000

SparseMatrix<double,RowMajor> *FGBGMatte::GetMatte() {
  // From "Scalable Matting: A Sub-linear Approach," by Lee and Wu.

  SparseMatrix<double,RowMajor> *L = source_->GetLaplacian();

  SparseMatrix<double,RowMajor> *g = guess_->GetGreyscaleVector();

  SparseMatrix<double,RowMajor> f(L->rows(), 1);
  SparseMatrix<double,RowMajor> C(L->rows(), L->cols());

  for (int i = 0; i < g->rows(); i++) {
    double guess = g->coeff(i,0);
    if (guess >= 1 || guess == 0) {
      f.insert(i,0) = guess;
      C.insert(i,i) = 1;
    } else {
      f.insert(i,0) = 0;
    }
  }

  // return g;

  f *= GAMMA;

  SparseMatrix<double,RowMajor> A(L->rows(), L->cols());
  A = (*L) + C * GAMMA;

  GridSolver solver(&A, &f);

  for (int i = 0; i < 100; i++) {
    std::cout << "Iter: " << i << ";    Resid. Norm: " << (A * *g).norm() << std::endl;
    solver.solve(g, 0); // For now, I won't use multigrid so
                        // that I can get the matting part working correctly on its own.
  }

  std::cout << "NORM: " << g->norm() << std::endl;

  // *g *= 1000;

  // for (int i = 0; i < g->rows(); i++) {
  //   if (g->coeff(i,0) < .001) {
  //     g->coeffRef(i,0) = 0;
  //   }
  // }

  SparseMatrix<double,RowMajor> *matte = new SparseMatrix<double,RowMajor>(source_->GetW(), source_->GetH());
  for (int i = 0; i < matte->rows(); ++i) {
    for (int j = 0; j < matte->cols(); ++j) {
      matte->insert(i,j) = g->coeff(i * source_->GetH() + j, 0);
    }
  }

  return matte;
}