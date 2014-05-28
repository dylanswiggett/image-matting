#include <iostream>
#include "GridSolver.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

#define SIZE 3000

int main(int argc, char** argv) {
  compressed_matrix<double> *m1;
  compressed_matrix<double> *m2;

  m1 = new compressed_matrix<double>(SIZE,SIZE,SIZE*SIZE);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      (*m1)(i,j) = i * SIZE + j;
    }
  }

  m2 = new compressed_matrix<double>(SIZE,SIZE,SIZE*SIZE);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      (*m2)(i,j) = i * SIZE + j + SIZE * SIZE;
    }
  }

  // std::cout << *m1 << std::endl;
  // std::cout << *m2 << std::endl;
  // std::cout << prod(*m1,*m2) << std::endl;
}