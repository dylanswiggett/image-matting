#include <iostream>
#include <string>
#include "GridSampler.hpp"
#include "ImageManager.hpp"
#include "FGBGMatte.hpp"
#include "SDL/SDL.h"

#include <eigen3/Eigen/Sparse>

using namespace Eigen;

int main(int argc, char** argv) {
  SDL_Init(SDL_INIT_EVERYTHING);

  std::string image_path = "test_img/toast.bmp";
  std::string img_guess_path = "test_img/toast_guess.bmp";
  std::string img_out_path = "test_img/toast_result.bmp";

  if (argc == 4) {
    image_path = argv[1];
    img_guess_path = argv[2];
    img_out_path = argv[3];
  }

  ImageManager img(image_path);
  ImageManager img_guess(img_guess_path);

  FGBGMatte matte(&img, &img_guess);

  ImageManager matte_img(matte.GetMatte());

  matte_img.SaveTo(img_out_path);
  return 0;
}
