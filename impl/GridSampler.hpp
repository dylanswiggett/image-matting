#ifndef _GRID_SAMPLER_HPP_
#define _GRID_SAMPLER_HPP_

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class GridSampler {
 public:
  GridSampler(const compressed_matrix<double>* mat) : mat_(mat) {};
  compressed_matrix<double>* upsample();
  compressed_matrix<double>* downsample();
 private:
  const compressed_matrix<double>* mat_;
};

#endif  // _GRID_SAMPLER_HPP_
