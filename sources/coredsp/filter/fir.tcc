#pragma once
#include "../filter.h"
#include "coreutil/restrict_ptr.h"
#include <algorithm>

namespace coredsp {

template <unsigned N, class V>
FIR<N, V>::FIR() {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 1) & ~(v - 1);

  c_.reset(new (coreutil::aligned(coreutil::simd_alignment)) real_type[nn]());
  x_.reset(new real_type[2 * nn]());
}

template <unsigned N, class V>
template <class Rc>
void FIR<N, V>::coefs(const Rc *coefs) {
  std::copy_n(coefs, N, c_.get());
}

template <unsigned N, class V>
inline void FIR<N, V>::reset() {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 1) & ~(v - 1);

  std::fill_n(x_.get(), 2 * nn, 0);
}

template <unsigned N, class V>
inline void FIR<N, V>::in(real_type x) {
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 1) & ~(v - 1);

  unsigned i = i_ = (i_ + nn - 1) % nn;
  x_[i] = x_[i + nn] = x;
}

template <unsigned N, class V>
inline auto FIR<N, V>::out() const -> real_type {
  // return impl_scalar();
  return impl_simd();
}

template <unsigned N, class V>
inline auto FIR<N, V>::tick(real_type x) -> real_type {
  in(x);
  return out();
}

template <unsigned N, class V>
inline auto FIR<N, V>::impl_scalar() const -> real_type {
  const unsigned i = i_;
  coreutil::restrict_ptr<const real_type> x = x_.get(), c = c_.get();

  real_type r = 0;
  for (unsigned j = 0; j < N; ++j)
    r += x[i + j] * c[j];
  return r;
}

template <unsigned N, class V>
inline auto FIR<N, V>::impl_simd() const -> real_type {
  const unsigned i = i_;
  constexpr unsigned v = coreutil::simd_size<V>;
  constexpr unsigned nn = (N + v - 1) & ~(v - 1);
  coreutil::restrict_ptr<const real_type> x = x_.get(), c = c_.get();

  vector_type rs {};
  for (unsigned j = 0; j < nn; j += v) {
    vector_type cs = *(const vector_type *)&c[j];
    vector_type xs = coreutil::simd_loadu((const vector_type *)&x[i + j]);
    rs = coreutil::simd_mul_add(xs, cs, rs);
  }
  return coreutil::simd_sum(rs);
}

}  // namespace coredsp
