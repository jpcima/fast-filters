#pragma once
#include "../filter.h"
#include "coreutil/restrict_ptr.h"
#include <algorithm>
#include <stdexcept>

namespace coredsp {

template <class V>
inline IIRg<V>::IIRg(unsigned n) : n_(n) {
  if (n < 2)
    throw std::runtime_error("invalid number of filter coefficients");

  constexpr unsigned v = coreutil::simd_size<V>;
  const unsigned nn = (n + v - 2) & ~(v - 1);

  b_.reset(new (coreutil::aligned(coreutil::simd_alignment)) real_type[nn]());
  a_.reset(new (coreutil::aligned(coreutil::simd_alignment)) real_type[nn]());
  x_.reset(new real_type[2 * nn]());
  y_.reset(new real_type[2 * nn]());
}

template <class V>
template <class Rc>
void IIRg<V>::coefs(const Rc *bcf, const Rc *acf) {
  const unsigned n = n_;
  b0_ = bcf[0];
  a0_ = acf[0];
  std::copy_n(bcf + 1, n - 1, b_.get());
  std::copy_n(acf + 1, n - 1, a_.get());
}

template <class V>
inline void IIRg<V>::reset() {
  constexpr unsigned v = coreutil::simd_size<V>;
  const unsigned nn = (n_ + v - 2) & ~(v - 1);

  std::fill_n(x_.get(), 2 * nn, 0);
  std::fill_n(y_.get(), 2 * nn, 0);
}

template <class V>
inline auto IIRg<V>::tick(real_type in) -> real_type {
  // return impl_scalar(in);
  return impl_simd(in);
}

template <class V>
inline auto IIRg<V>::impl_scalar(real_type in) -> real_type {
  constexpr unsigned v = coreutil::simd_size<V>;
  const unsigned n = n_;
  const unsigned nn = (n + v - 2) & ~(v - 1);
  coreutil::restrict_ptr<real_type> x = x_.get(), y = y_.get();
  coreutil::restrict_ptr<const real_type> bcf = b_.get(), acf = a_.get();

  real_type r = in * b0_;
  unsigned i = (i_ + nn - 1) % nn;
  x[i] = x[i + nn] = in;

  for (unsigned j = 0; j < n - 1; ++j)
    r += x[i + j + 1] * bcf[j] - y[i + j + 1] * acf[j];

  r /= a0_;
  y[i] = y[i + nn] = r;
  i_ = i;
  return r;
}

template <class V>
inline auto IIRg<V>::impl_simd(real_type in) -> real_type {
  constexpr unsigned v = coreutil::simd_size<V>;
  const unsigned n = n_;
  const unsigned nn = (n + v - 2) & ~(v - 1);
  coreutil::restrict_ptr<real_type> x = x_.get(), y = y_.get();
  coreutil::restrict_ptr<const real_type> bcf = b_.get(), acf = a_.get();

  real_type r = in * b0_;
  unsigned i = (i_ + nn - 1) % nn;
  x[i] = x[i + nn] = in;

  vector_type rs {};
  for (unsigned j = 0; j < nn; j += v) {
    vector_type bs = *(const vector_type *)&bcf[j];
    vector_type xs = coreutil::simd_loadu((const vector_type *)&x[i + j + 1]);
    vector_type as = *(const vector_type *)&acf[j];
    vector_type ys = coreutil::simd_loadu((const vector_type *)&y[i + j + 1]);
    rs = coreutil::simd_mul_add(xs, bs, rs);
    rs = coreutil::simd_mul_neg_add(ys, as, rs);
  }
  r += coreutil::simd_sum(rs);

  r /= a0_;
  y[i] = y[i + nn] = r;
  i_ = i;
  return r;
}

}  // namespace coredsp
