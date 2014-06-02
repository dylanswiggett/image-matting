#include <iostream>
#include "GridSolver.hpp"
#include "GridSampler.hpp"
#include "ImageManager.hpp"
#include "FGBGMatte.hpp"
#include "SDL/SDL.h"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

// #define SIZE 10000

int main(int argc, char** argv) {
  SDL_Init(SDL_INIT_EVERYTHING);
  // SparseMatrix<double,RowMajor> *m;

  // m = new SparseMatrix<double,RowMajor>(SIZE,SIZE);
  // m->reserve(SIZE*3);
  // for (int i = 0; i < SIZE; i++) {
  //   for (int j = 0; j < SIZE; j++) {
  //     if (i == j)
  //       m->insert(i,j) = 2;
  //     else if (i - j == 1 || i - j == -1)
  //       m->insert(i,j) = 1;
  //   }
  // }
  // m->makeCompressed();

  ImageManager img("test_img/pooh_med.bmp");
  ImageManager img_guess("test_img/pooh_med_guess.bmp");

  ImageManager(GridSampler(img_guess.GetGreyscaleMatrix()).upsample(img_guess.GetW() * 2, img_guess.GetH() * 2)).SaveTo("test_img/pooh_med_guess_copy.bmp");

  // ImageManager(img.GetLaplacian()).SaveTo("test_img/pooh_med_laplacian.bmp");

  FGBGMatte matte(&img, &img_guess);

  ImageManager matte_img(matte.GetMatte());

  matte_img.SaveTo("test_img/pooh_med_result.bmp");
  // ImageManager(img.GetGreyscaleMatrix()).SaveTo("test_img/small_2color_resave.bmp");
  // ImageManager(img_guess.GetGreyscaleMatrix()).SaveTo("test_img/small_2color_guess_resave.bmp");

  // std::cout << *(matte.GetMatte()) << std::endl;

  // std::cout << *(img.GetLaplacian()) << std::endl;

  // std::cout << *(GridSampler(GridSampler(m).upsample(7,7)).downsample()) << std::endl;

  return 0;

  // std::cout << "'A' Matrix:\n--------------------------------\n" << *m << "\n--------------------------------" << std::endl;

  // SparseMatrix<double,RowMajor> *expected;

  // expected = new SparseMatrix<double,RowMajor>(SIZE,1);

  // for (int i = 0; i < SIZE; i++)
  //   (*expected).insert(i,0) = 0;

  // GridSolver solver(m,expected);

  // SparseMatrix<double,RowMajor> *guess = new SparseMatrix<double,RowMajor>(SIZE,1);

  // for (int i = 0; i < SIZE; i++)
  //   (*guess).insert(i,0) = 1;

  // std::cout << "Stepping once... " << std::endl;

  // for (int i = 0; i < 4; i++)
  //   solver.solve(guess, 1);

  // std::cout << "Done!\nFound solution:\n" << *guess << "\nWhich has error:" << std::endl;

  // std::cout << (*m) * (*guess) - *expected << std::endl;

  // delete m;
  // delete expected;
  // delete guess;
}
