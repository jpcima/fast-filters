# fast-filters
Implementation of FIR and IIR filters optimized for SIMD processing

## Benchmarks

The following benchmarks have been computed on x86_64 Linux with Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz.

| AVX implementation                                   | SSE implementation                                   |
| ---------------------------------------------------- | ---------------------------------------------------- |
| ![fir-avx-float](benchmark-data/fir-avx-float.png)   | ![fir-sse-float](benchmark-data/fir-sse-float.png)   |
| ![fir-avx-double](benchmark-data/fir-avx-double.png) | ![fir-sse-double](benchmark-data/fir-sse-double.png) |
| ![iir-avx-float](benchmark-data/iir-avx-float.png)   | ![iir-sse-float](benchmark-data/iir-sse-float.png)   |
| ![iir-avx-double](benchmark-data/iir-avx-double.png) | ![iir-sse-double](benchmark-data/iir-sse-double.png) |
