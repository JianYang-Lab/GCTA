# GCTA
GCTA (Genome-wide Complex Trait Analysis) is a software package, which was initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait but has been extensively extended for many other analyses of data from genome-wide association studies (GWASs). Please see the software website through the link below for more information.

Software website: https://yanglab.westlake.edu.cn/software/gcta/
License: GPLv3 (some parts of the code are released under LGPL as detailed in the files).


## Credits  
Jian Yang developed the original version (before v1.90) of the software (with supports from Peter Visscher, Mike Goddard and Hong Lee) and currently maintains the software.

Zhili Zheng programmed the fastGWA, fastGWA-GLMM and fastGWA-BB modules, rewrote the I/O and GRM modules, improved the GREML and bivariate GREML modules, extended the PCA module, and improved the SBLUP module.  

Zhihong Zhu programmed the mtCOJO and GSMR modules and improved the COJO module.  

Longda Jiang and Hailing Fang developed the ACAT-V module.  

Jian Zeng rewrote the GCTA-HEreg module.  

Andrew Bakshi contributed to the GCTA-fastBAT module.

Angli Xue improved the GSMR module.

Robert Maier improved the GCTA-SBLUP module.

Wujuan Zhong and Judong Shen programmed the fastGWA-GE module. 

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.


## Compilation

#### Requirements

1. Currently only x86\_64-based operating systems are supported.
2. [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) 2017 or above (only needed when building on x86\-64 machines)
3. OpenBLAS (only needed when building on AArch64 machines)
4. Eigen == 3.3.7 (there are bugs in the new version of Eigen)
5. CMake >= 3.1
6. BOOST >= 1.4
7. zlib >= 1.2.11 (old zlib may cause an error of bgen file decompression)
8. sqlite3 >= 3.31.1
9. zstd >= 1.4.4
10. [Spectra](https://spectralib.org/) >= 0.8.1
11. gsl (GNU scientific library)

Most of the dependencies above will be downloaded by CMake automatically. You only need to install the `gsl` and `Intel MKL` manually.

#### Linux

1. Kernel version >= 2.6.28 (otherwise the Intel MKL library doesn't work).
2. GCC version >= 6.1 with C++ 11 support.

#### Before compilation 

Update [plink_ng](https://github.com/chrchang/plink-ng) submodule first.

```sh
git submodule update --init
```

On Windows, apply the patch under the `third_party` directory to the `plink-ng`.

#### Build

##### CMake Configuration

On MacOS and Linux, use following command to generate the build system:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja -B build/Release -S .
```

On Windows, you should use the toolchain file in `cmake/win-toolchain.cmake`:

``` sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="cmake/win-toolchain.cmake" -G Ninja -B build/Release -S .
```

##### Compile

```sh
cmake --build build/Release
```

The executable binary will be generated under `build/Release`.

