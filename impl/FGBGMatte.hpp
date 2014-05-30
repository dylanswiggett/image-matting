#ifndef _FGBGMATTE_HPP_
#define _FGBGMATTE_HPP_

#include "./ImageManager.hpp"

class FGBGMatte {
 public:
  FGBGMatte(ImageManager source, ImageManager guess);
  ImageManager GetMatte();
 private:
  ImageManager source_ guess_;
};

#endif  // _FGBGMATTE_HPP_