set(CMAKE_SYSTEM_NAME Linux)

# module load gcc/8.4.0 intelmkl gsl
set(CMAKE_C_COMPILER /soft/compiler/gcc/gcc-8.4.0/bin/gcc)
set(CMAKE_CXX_COMPILER /soft/compiler/gcc/gcc-8.4.0/bin/g++)
set(CMAKE_PREFIX_PATH "/soft/mathlib/gsl-2.6;")
set(MKL_DIR /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/cmake/mkl)

