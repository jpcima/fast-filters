#include "coredsp/filter.h"
#include "coredsp/noise.h"
#include "coreutil/cpu.h"
#include <boost/preprocessor.hpp>
#include <functional>
#include <chrono>
#include <vector>

typedef std::chrono::steady_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;

static unsigned iteration_count = 64 * 1024;

struct Benchmark {
  const char *name {};
  std::function<float()> fn;
  double run() const;
};

double Benchmark::run() const {
  TimePoint tstart = Clock::now();
  this->fn();
  Duration dur = Clock::now() - tstart;
  auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(dur);
  return ns.count() * 1e-9;
}

#define IIR_NCOEF_MIN 2
#define IIR_NCOEF_LIMIT 256

static constexpr unsigned niir = IIR_NCOEF_LIMIT - IIR_NCOEF_MIN;

#define DECL_IIR_FLT(z, n, t)                               \
  static coredsp::IIR<n, coreutil::simd_t<float>> iirf##n;  \
  static coredsp::IIR<n, coreutil::simd_t<double>> iird##n;
#define DECL_IIRG_FLT(z, n, t)                                 \
  static coredsp::IIRg<coreutil::simd_t<float>> iirgf##n(n);   \
  static coredsp::IIRg<coreutil::simd_t<double>> iirgd##n(n);
BOOST_PP_REPEAT_FROM_TO(IIR_NCOEF_MIN, IIR_NCOEF_LIMIT, DECL_IIR_FLT,);
BOOST_PP_REPEAT_FROM_TO(IIR_NCOEF_MIN, IIR_NCOEF_LIMIT, DECL_IIRG_FLT,);

static Benchmark biirf_simd[niir], biirf_scalar[niir];
static Benchmark biirgf_simd[niir], biirgf_scalar[niir];
static Benchmark biird_simd[niir], biird_scalar[niir];
static Benchmark biirgd_simd[niir], biirgd_scalar[niir];

template <class T> std::function<float()> create_iir_simd_fn(T &filter) {
  return [&filter]() -> float {
    float out {};
    coredsp::WhiteNoise source;
    for (unsigned i = 0; i < iteration_count; ++i)
      out = filter.impl_simd(source.tick());
    return out;
  };
}

template <class T> std::function<float()> create_iir_scalar_fn(T &filter) {
  return [&filter]() -> float {
    float out {};
    coredsp::WhiteNoise source;
    for (unsigned i = 0; i < iteration_count; ++i)
      out = filter.impl_scalar(source.tick());
    return out;
  };
}

int main() {
  coreutil::disable_denormals();

#define DEF_IIR_BENCHMARK(z, n, t)                                      \
  biirf_simd[n - IIR_NCOEF_MIN] = Benchmark{"IIR<" #n "> simd float", create_iir_simd_fn(iirf##n)}; \
  biirf_scalar[n - IIR_NCOEF_MIN] = Benchmark{"IIR<" #n "> scalar float", create_iir_scalar_fn(iirf##n)}; \
  biirgf_simd[n - IIR_NCOEF_MIN] = Benchmark{"IIRg<" #n "> simd float", create_iir_simd_fn(iirgf##n)}; \
  biirgf_scalar[n - IIR_NCOEF_MIN] = Benchmark{"IIRg<" #n "> scalar float", create_iir_scalar_fn(iirgf##n)}; \
  biird_simd[n - IIR_NCOEF_MIN] = Benchmark{"IIR<" #n "> simd double", create_iir_simd_fn(iird##n)}; \
  biird_scalar[n - IIR_NCOEF_MIN] = Benchmark{"IIR<" #n "> scalar double", create_iir_scalar_fn(iird##n)}; \
  biirgd_simd[n - IIR_NCOEF_MIN] = Benchmark{"IIRg<" #n "> simd double", create_iir_simd_fn(iirgd##n)}; \
  biirgd_scalar[n - IIR_NCOEF_MIN] = Benchmark{"IIRg<" #n "> scalar double", create_iir_scalar_fn(iirgd##n)};

  BOOST_PP_REPEAT_FROM_TO(IIR_NCOEF_MIN, IIR_NCOEF_LIMIT, DEF_IIR_BENCHMARK,);

  setlinebuf(stdout);

  for (unsigned i = 0; i < niir; ++i) {
    double tsimdf = biirf_simd[i].run();
    double tscalarf = biirf_scalar[i].run();
    double tgsimdf = biirgf_simd[i].run();
    double tgscalarf = biirgf_scalar[i].run();
    double tsimdd = biird_simd[i].run();
    double tscalard = biird_scalar[i].run();
    double tgsimdd = biirgd_simd[i].run();
    double tgscalard = biirgd_scalar[i].run();

    printf("%u %f %f %f %f %f %f %f %f\n", i + IIR_NCOEF_MIN,
           tsimdf, tscalarf, tgsimdf, tgscalarf,
           tsimdd, tscalard, tgsimdd, tgscalard);
  }

  return 0;
}
