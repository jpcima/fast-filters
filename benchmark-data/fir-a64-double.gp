set title "FIR double precision A64 - 64k iterations"

set xlabel "Length of filter"
set ylabel "Time (s)"

plot "fir-a64.dat" u 1:6 w lines t "SIMD compile-time N", \
     "" u 1:7 w lines t "Scalar compile-time N", \
     "" u 1:8 w lines t "SIMD run-time N", \
     "" u 1:9 w lines t "Scalar run-time N"
