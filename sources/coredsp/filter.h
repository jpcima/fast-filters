#pragma once
#include "coreutil/simd.h"
#include "coreutil/memory.h"
#include <memory>

namespace coredsp {

//------------------------------------------------------------------------------
// FIR filter with compile-time length N
//------------------------------------------------------------------------------
template <unsigned N, class V = coreutil::simd_t<float>>
struct FIR {
  static_assert(N >= 2, "invalid number of filter coefficients");

  typedef V vector_type;
  typedef coreutil::simd_element_type<V> real_type;

  FIR();
  template <class Rc> void coefs(const Rc *coefs);
  void reset();
  void in(real_type x);
  real_type out() const;
  real_type tick(real_type x);

  FIR(FIR &&) = default;
  FIR &operator=(FIR &&) = default;

 private:
  unsigned i_ = {};
  std::unique_ptr<real_type[]> x_;
  std::unique_ptr<real_type[], coreutil::aligned_delete<real_type>> c_;

 public:
  real_type impl_scalar() const;
  real_type impl_simd() const;
};

//------------------------------------------------------------------------------
// FIR filter with run-time length
//------------------------------------------------------------------------------
template <class V = coreutil::simd_t<float>>
struct FIRg {
  typedef V vector_type;
  typedef coreutil::simd_element_type<V> real_type;

  explicit FIRg(unsigned n);
  template <class Rc> void coefs(const Rc *coefs);
  void reset();
  void in(real_type x);
  real_type out() const;
  real_type tick(real_type x);

  FIRg(FIRg &&) = default;
  FIRg &operator=(FIRg &&) = default;

 private:
  const unsigned n_ {};
  unsigned i_ {};
  std::unique_ptr<real_type[]> x_;
  std::unique_ptr<real_type[], coreutil::aligned_delete<real_type>> c_;

 public:
  real_type impl_scalar() const;
  real_type impl_simd() const;
};

//------------------------------------------------------------------------------
// IIR filter with compile-time length N
//------------------------------------------------------------------------------
template <unsigned N, class V = coreutil::simd_t<double>>
struct IIR {
  static_assert(N >= 2, "invalid number of filter coefficients");

  typedef V vector_type;
  typedef coreutil::simd_element_type<V> real_type;

  IIR();
  template <class Rc> void coefs(const Rc *bcf, const Rc *acf);
  void reset();
  real_type tick(real_type in);

  IIR(IIR &&) = default;
  IIR &operator=(IIR &&) = default;

 private:
  unsigned i_ = {};
  std::unique_ptr<real_type[]> x_, y_;
  real_type b0_ = 1, a0_ = 1;
  std::unique_ptr<real_type[], coreutil::aligned_delete<real_type>> b_, a_;

 public:
  real_type impl_scalar(real_type in);
  real_type impl_simd(real_type in);
};

//------------------------------------------------------------------------------
// FIR filter with run-time length
//------------------------------------------------------------------------------
template <class V = coreutil::simd_t<double>>
struct IIRg {
  typedef V vector_type;
  typedef coreutil::simd_element_type<V> real_type;

  explicit IIRg(unsigned n);
  template <class Rc> void coefs(const Rc *bcf, const Rc *acf);
  void reset();
  real_type tick(real_type in);

  IIRg(IIRg &&) = default;
  IIRg &operator=(IIRg &&) = default;

 private:
  const unsigned n_ {};
  unsigned i_ = {};
  std::unique_ptr<real_type[]> x_, y_;
  real_type b0_ = 1, a0_ = 1;
  std::unique_ptr<real_type[], coreutil::aligned_delete<real_type>> b_, a_;

 public:
  real_type impl_scalar(real_type in);
  real_type impl_simd(real_type in);
};

}  // namespace coredsp

#include "filter/fir.tcc"
#include "filter/firg.tcc"
#include "filter/iir.tcc"
#include "filter/iirg.tcc"
