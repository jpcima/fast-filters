#pragma once
#include "../filter.h"
#include "coreutil/restrict_ptr.h"
#include <algorithm>

namespace coredsp {

template <unsigned N, class V>
IIR<N, V>::IIR() {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 2) & ~(v - 1);

  b_.reset(new (coreutil::aligned(coreutil::simd_alignment)) real_type[nn]());
  a_.reset(new (coreutil::aligned(coreutil::simd_alignment)) real_type[nn]());
  x_.reset(new real_type[2 * nn]());
  y_.reset(new real_type[2 * nn]());
}

template <unsigned N, class V>
template <class Rc>
void IIR<N, V>::coefs(const Rc *bcf, const Rc *acf) {
  b0_ = bcf[0];
  a0_ = acf[0];
  std::copy_n(bcf + 1, N - 1, b_.get());
  std::copy_n(acf + 1, N - 1, a_.get());
}

template <unsigned N, class V>
inline void IIR<N, V>::reset() {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 2) & ~(v - 1);

  std::fill_n(x_, 2 * nn, 0);
  std::fill_n(y_, 2 * nn, 0);
}

template <unsigned N, class V>
inline auto IIR<N, V>::tick(real_type in) -> real_type {
  // return impl_scalar(in);
  return impl_simd(in);
}

template <unsigned N, class V>
inline auto IIR<N, V>::impl_scalar(real_type in) -> real_type {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 2) & ~(v - 1);
  coreutil::restrict_ptr<real_type> x = x_.get(), y = y_.get();
  coreutil::restrict_ptr<const real_type> bcf = b_.get(), acf = a_.get();

  real_type r = in * b0_;
  unsigned i = (i_ + nn - 1) % nn;
  x[i] = x[i + nn] = in;

  for (unsigned j = 0; j < N - 1; ++j)
    r += x[i + j + 1] * bcf[j] - y[i + j + 1] * acf[j];

  r /= a0_;
  y[i] = y[i + nn] = r;
  i_ = i;
  return r;
}

template <unsigned N, class V>
inline auto IIR<N, V>::impl_simd(real_type in) -> real_type {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 2) & ~(v - 1);
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
