# QC-MDPC decoder

This is an implementation of a few QC-MDPC decoding algorithms.

See <https://tel.archives-ouvertes.fr/tel-03254461/>.


## Usage

```sh
./qcmdpc_decoder [OPTIONS]

-i, --max-iter         maximum number of iterations
-N, --rounds           number of rounds to perform
-T, --threads          number of threads to use
-q, --quiet            do not regularly output results (only on SIGHUP)
```

It generates QC-MDPC decoding instances then tries to decode them. For each
instance, a random parity check matrix and a random error vector are generated
then the corresponding syndrome is computed.

Every 5 seconds, it prints the number of instances generated and the
distribution of the number of iterations it took to decode.

Unless a number of rounds is specified, it will only stop on SIGINT (Ctrl+C) or
SIGTERM.


## Example

```sh
$ ./qcmdpc_decoder -T4 -N100000
-DINDEX=2 -DBLOCK_LENGTH=10007 -DBLOCK_WEIGHT=71 -DERROR_WEIGHT=134 -DOUROBOROS=0 -DWEAK=0 -DWEAK_P=0 -DERROR_FLOOR=0 -DERROR_FLOOR_P=0 -DTHRESHOLD_C0=13.53 -DTHRESHOLD_C1=0.0069722 -DALGO=GRAY_BGF
58907 5:20339 7:38396 9:163 11:6 >100:4
100000 5:34466 7:65226 9:294 11:10 >100:4
```

Out of 100000 instances, 34466 were decoded in 5 iterations, 65226 in 7
iterations, 294 in 9 iterations, 10 in 11 iterations and 4 failed to decode
in at most 100 iterations.


## Parameters

For speed, parameters are chosen at compile time. They are:
- `INDEX`: number of circulant blocks in the parity check matrix,
- `BLOCK_LENGTH`,
- `BLOCK_WEIGHT`: column weight of the parity check matrix,
- `ERROR_WEIGHT`,
- `OUROBOROS` (0 or 1): noisy syndrome variant,
- `WEAK` (0-3 depending on the type): weak key generation,
- `WEAK_P`: number of successive ones for Type I, maximum multiplicity in the distance spectrum for types II and III,
- `ERROR_FLOOR` (0-3) error patterns close to (1) (d, d) near-codewords, (2) (2d, ~2d) near-codewords, (3) codewords,
- `ERROR_FLOOR_P`: number of intersections between error patterns and (near-)codewords.

Algorithm and their respective parameters can be chosen among:
- `ALGO = BACKFLIP`: Backflip with an affine ttl
    * `TTL_C0`, `TTL_C1`: affine ttl function coefficients
    * `TTL_SATURATE`: ttl saturation value
- `ALGO = BACKFLIP2`: Backflip with multiple thresholds
    * `THRESHOLD_A0`, `THRESHOLD_A1`, `THRESHOLD_A2`, `THRESHOLD_A3`, `THRESHOLD_A4`: ttl alpha values
    * `TTL_SATURATE`: ttl saturation value
- `ALGO = BP`: belief propagation
    * `BP_SCALE`: scaling factor in messages from check to variable nodes
    * `BP_SATURATE` messages saturation value
- `ALGO = CLASSIC`: classic bit-flipping algorithm
- `ALGO = GRAY_B |  GRAY_BGF | GRAY_BGB | GRAY_BG`
    * `THRESHOLD_C0`, `THRESHOLD_C1`: affine threshold function coefficients
- `ALGO = SBS`: step-by-step algorithm
- `ALGO = SORT`: sorted gray algorithm
    * `GRAY_SIZE`

Build with, for example:
```sh
$ cmake -B build/ -DINDEX=2 -DBLOCK_LENGTH=12323 -DBLOCK_WEIGHT=71 -DERROR_WEIGHT=134 -DOUROBOROS=0 -DWEAK=0 -DWEAK_P=0 -DERROR_FLOOR=0 -DERROR_FLOOR_P=0 -DTHRESHOLD_C0=13.53 -DTHRESHOLD_C1=0.0069722 -DALGO=GRAY_BGB && cmake --build build/
```

Executable name is `qcmdpc_decoder` located in the `build` directory.


## Parameters presets

Presets are available for BIKE parameters:
- `PRESET_CPA = 128`,
- `PRESET_CPA = 192`,
- `PRESET_CPA = 256`,
- `PRESET_CCA = 128`,
- `PRESET_CCA = 192`,
- `PRESET_CCA = 256`.

Build with, for example:
```sh
$ cmake -B build/ -DPRESET_CCA=192 -DALGO=BP && cmake --build build/
```


## Profile Guided Optimization

GCC and Clang do a good job at Profile Guided Optimization.
To use it, first compile with, for example:
```sh
$ cmake -B build/ -DINDEX=2 -DBLOCK_LENGTH=12323 -DBLOCK_WEIGHT=71 -DERROR_WEIGHT=134 -DALGO=BACKFLIP2 -DPGO=GEN && cmake --build build/
```

Run the program on a sample with, for example (for 8 iterations, 8 threads and
a sample of size 100000):
```sh
$ build/qcmdpc_decoder -i8 -T8 -N100000
```

Recompile to use PGO:
```sh
$ cmake -B build/ -DINDEX=2 -DBLOCK_LENGTH=12323 -DBLOCK_WEIGHT=71 -DERROR_WEIGHT=134 -DALGO=BACKFLIP2 -DPGO=USE && cmake --build build/
```


## AVX2

By default, the Makefile compiles the AVX2 version, if you do not have such an
instruction set, set the `AVX` option to `OFF`.
```sh
$ cmake -B build/ -DPRESET_CPA=256 -DALGO=CLASSIC -DAVX=OFF && cmake --build build/
```


# Scripts

A few scripts are included in the 'scripts/' directory.


## Summary of simulation

This script simply prints a summary of simulation data in a nice format. The
summary also includes a density estimate in the case of weak keys or error
floor patterns. For each number of iterations performed, the binary logarithm
of the DFR is given with Clopper-Pearson confidence intervals.

### Example of usage and output
```sh
$ cmake -B build/ -DBLOCK_LENGTH=9901 -DALGO=BACKFLIP2 && cmake --build build/
$ build/qcmdpc_decoder -T4 -N3000000 > 9901
$ python3 scripts/summary.py 9901 0.01
index        : 2
block_length : 9901
block_weight : 71
error_weight : 134
algo         : BACKFLIP2
threshold_a0 : 4
threshold_a1 : 1
threshold_a2 : 0.25
threshold_a3 : 0.0625
threshold_a4 : 0.015625
ttl_saturate : 5
dfr (with CI):
               3: -0.120 -0.119 -0.119
               4: -3.293 -3.286 -3.280
               5: -7.173 -7.148 -7.122
               6: -9.818 -9.754 -9.691
               7: -11.351 -11.244 -11.139
               8: -12.410 -12.257 -12.108
               9: -13.143 -12.947 -12.758
              10: -13.557 -13.332 -13.116
              11: -13.839 -13.592 -13.356
              12: -13.996 -13.735 -13.488
              13: -14.107 -13.837 -13.581
              14: -14.228 -13.947 -13.681
              15: -14.314 -14.025 -13.751
              16: -14.396 -14.099 -13.818
              17: -14.453 -14.150 -13.865
              18: -14.482 -14.177 -13.889
              19: -14.502 -14.195 -13.905
              20: -14.522 -14.213 -13.921
              21: -14.553 -14.240 -13.946
              22: -14.574 -14.259 -13.963
              23: -14.595 -14.278 -13.980
              25: -14.606 -14.288 -13.989
              30: -14.627 -14.307 -14.006
```

The Python 3 script requires gmpy2 and statsmodels.


## Extrapolation

This script returns the binary logarithm of the extrapolated DFR from two sets
of simulation data. It also includes a density estimate in the case of weak
keys or error floor patterns.

### Example of usage and output
```sh
$ cmake -B build/ -DBLOCK_LENGTH=9901 -DALGO=GRAY_BGF -DWEAK=1 -DWEAK_P=20 && cmake --build build/
[...]
$ build/qcmdpc_decoder -T4 -N100000 > 9901
$ cmake -B build/ -DBLOCK_LENGTH=10103 -DALGO=GRAY_BGF -DWEAK=1 -DWEAK_P=20 && cmake --build build/
[...]
$ build/qcmdpc_decoder -T4 -N100000 > 10103

$ python3 scripts/extrapolate.py 9901 10103 12323 0.01
Point 1
index        : 2
block_length : 9901
block_weight : 71
error_weight : 134
algo         : GRAY_BGF
threshold_c0 : 13.53
threshold_c1 : 0.0069722
weak         : 1
weak_p       : 20
density      : -120.170446669296

Point 2
index        : 2
block_length : 10103
block_weight : 71
error_weight : 134
algo         : GRAY_BGF
threshold_c0 : 13.53
threshold_c1 : 0.0069722
weak         : 1
weak_p       : 20
density      : -120.69547486997696

Extrapolated values:
block_length : 12323
density      : -125.85859847524915
dfr (with CI):
               4: -0.669 -0.666 -0.654
               5: -34.491 -34.295 -34.099
               6: -76.217 -75.496 -74.776
               7: -100.566 -98.586 -96.630
               8: -111.571 -108.057 -104.624
               9: -115.959 -111.363 -106.905
              10: -122.671 -116.795 -111.149
              11: -124.303 -117.884 -111.740
              12: -128.203 -120.983 -114.112
              13: -132.061 -124.019 -116.411
              14: -131.798 -123.756 -116.146
              15: -133.686 -125.254 -117.297
              16: -134.298 -125.723 -117.637
              17: -134.198 -125.621 -117.536
              18: -134.922 -126.194 -117.975
              19: -135.676 -126.789 -118.430
dfr+density (with CI):
               4: -126.528 -126.524 -126.512
               5: -160.349 -160.154 -159.957
               6: -202.076 -201.354 -200.634
               7: -226.425 -224.445 -222.489
               8: -237.429 -233.916 -230.483
               9: -241.818 -237.221 -232.764
              10: -248.530 -242.654 -237.008
              11: -250.162 -243.743 -237.599
              12: -254.061 -246.841 -239.970
              13: -257.920 -249.878 -242.270
              14: -257.656 -249.614 -242.004
              15: -259.545 -251.113 -243.156
              16: -260.157 -251.582 -243.496
              17: -260.057 -251.480 -243.394
              18: -260.781 -252.053 -243.834
              19: -261.535 -252.648 -244.288
```

The Python 3 script requires gmpy2 and scipy.


## Basic confidence intervals extrapolation

Confidence intervals for extrapolations can also be calculated "by hand".

```sh
$ python3 scripts/ci.py 12323 10037 66391 3747161784 10253 5 1445221866 0.01
-164.21 -146.20 -130.31
```

```sh
$ python3 scripts/ci.py 12323 10181 394 14576092619 10253 111 34283154045 0.01
-128.13 -116.22 -104.57
```

The Python 3 script requires gmpy2 and scipy.


# License

MIT (see file headers)
