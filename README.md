# fast-filters
Implementation of FIR and IIR filters optimized for SIMD processing

## Benchmarks

The following benchmarks have been computed on x86_64 Linux with Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz.

| AVX single precision                               | AVX double precision                                 |
| -------------------------------------------------- | ---------------------------------------------------- |
| ![fir-avx-float](benchmark-data/fir-avx-float.png) | ![fir-avx-double](benchmark-data/fir-avx-double.png) |
| ![iir-avx-float](benchmark-data/iir-avx-float.png) | ![iir-avx-double](benchmark-data/iir-avx-double.png) |

| SSE single precision                               | SSE double precision                                 |
| -------------------------------------------------- | ---------------------------------------------------- |
| ![fir-sse-float](benchmark-data/fir-sse-float.png) | ![fir-sse-double](benchmark-data/fir-sse-double.png) |
| ![iir-sse-float](benchmark-data/iir-sse-float.png) | ![iir-sse-double](benchmark-data/iir-sse-double.png) |

The following benchmarks have been computed on aarch64 Linux with Odroid-C2.

| A64 single precision                               | A64 double precision                                 |
| -------------------------------------------------- | ---------------------------------------------------- |
| ![fir-a64-float](benchmark-data/fir-a64-float.png) | ![fir-a64-double](benchmark-data/fir-a64-double.png) |
| ![iir-a64-float](benchmark-data/iir-a64-float.png) | ![iir-a64-double](benchmark-data/iir-a64-double.png) |
