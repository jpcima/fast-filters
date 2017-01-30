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

#define FIR_NCOEF_MIN 2
#define FIR_NCOEF_LIMIT 256

static constexpr unsigned nfir = FIR_NCOEF_LIMIT - FIR_NCOEF_MIN;

#define DECL_FIR_FLT(z, n, t)                               \
  static coredsp::FIR<n, coreutil::simd_t<float>> firf##n;  \
  static coredsp::FIR<n, coreutil::simd_t<double>> fird##n;
#define DECL_FIRG_FLT(z, n, t)                                 \
  static coredsp::FIRg<coreutil::simd_t<float>> firgf##n(n);   \
  static coredsp::FIRg<coreutil::simd_t<double>> firgd##n(n);
BOOST_PP_REPEAT_FROM_TO(FIR_NCOEF_MIN, FIR_NCOEF_LIMIT, DECL_FIR_FLT,);
BOOST_PP_REPEAT_FROM_TO(FIR_NCOEF_MIN, FIR_NCOEF_LIMIT, DECL_FIRG_FLT,);

static Benchmark bfirf_simd[nfir], bfirf_scalar[nfir];
static Benchmark bfirgf_simd[nfir], bfirgf_scalar[nfir];
static Benchmark bfird_simd[nfir], bfird_scalar[nfir];
static Benchmark bfirgd_simd[nfir], bfirgd_scalar[nfir];

template <class T> std::function<float()> create_fir_simd_fn(T &filter) {
  return [&filter]() -> float {
    float out {};
    coredsp::WhiteNoise source;
    for (unsigned i = 0; i < iteration_count; ++i) {
      filter.in(source.tick());
      out = filter.impl_simd();
    }
    return out;
  };
}

template <class T> std::function<float()> create_fir_scalar_fn(T &filter) {
  return [&filter]() -> float {
    float out {};
    coredsp::WhiteNoise source;
    for (unsigned i = 0; i < iteration_count; ++i) {
      filter.in(source.tick());
      out = filter.impl_scalar();
    }
    return out;
  };
}

int main() {
  coreutil::disable_denormals();

#define DEF_FIR_BENCHMARK(z, n, t)                                      \
  bfirf_simd[n - FIR_NCOEF_MIN] = Benchmark{"FIR<" #n "> simd float", create_fir_simd_fn(firf##n)}; \
  bfirf_scalar[n - FIR_NCOEF_MIN] = Benchmark{"FIR<" #n "> scalar float", create_fir_scalar_fn(firf##n)}; \
  bfirgf_simd[n - FIR_NCOEF_MIN] = Benchmark{"FIRg<" #n "> simd float", create_fir_simd_fn(firgf##n)}; \
  bfirgf_scalar[n - FIR_NCOEF_MIN] = Benchmark{"FIRg<" #n "> scalar float", create_fir_scalar_fn(firgf##n)}; \
  bfird_simd[n - FIR_NCOEF_MIN] = Benchmark{"FIR<" #n "> simd double", create_fir_simd_fn(fird##n)}; \
  bfird_scalar[n - FIR_NCOEF_MIN] = Benchmark{"FIR<" #n "> scalar double", create_fir_scalar_fn(fird##n)}; \
  bfirgd_simd[n - FIR_NCOEF_MIN] = Benchmark{"FIRg<" #n "> simd double", create_fir_simd_fn(firgd##n)}; \
  bfirgd_scalar[n - FIR_NCOEF_MIN] = Benchmark{"FIRg<" #n "> scalar double", create_fir_scalar_fn(firgd##n)};

  BOOST_PP_REPEAT_FROM_TO(FIR_NCOEF_MIN, FIR_NCOEF_LIMIT, DEF_FIR_BENCHMARK,);

  setlinebuf(stdout);

  for (unsigned i = 0; i < nfir; ++i) {
    double tsimdf = bfirf_simd[i].run();
    double tscalarf = bfirf_scalar[i].run();
    double tgsimdf = bfirgf_simd[i].run();
    double tgscalarf = bfirgf_scalar[i].run();
    double tsimdd = bfird_simd[i].run();
    double tscalard = bfird_scalar[i].run();
    double tgsimdd = bfirgd_simd[i].run();
    double tgscalard = bfirgd_scalar[i].run();

    printf("%u %f %f %f %f %f %f %f %f\n", i + FIR_NCOEF_MIN,
           tsimdf, tscalarf, tgsimdf, tgscalarf,
           tsimdd, tscalard, tgsimdd, tgscalard);
  }

  return 0;
}
