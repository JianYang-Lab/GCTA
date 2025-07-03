# Build GCTA

This document described general procedure to build GCTA, for more detail please read our document for developer.

## Linux

1. Requirements

    * libs

        1. intel MKL or OpenBLAS(until now only aarch64 using OpenBLAS and other using MKL)

        2. Eigen == 3.3.7

        3. Spectra >= 0.8.1

        4. Boost >= 1.4

        5. zlib >= 1.2.11

        6. zstd >= 1.4.4

        7. gsl(only version 2.7 was tested)

        8. sqlite3 >= 3.31.1

        You can istall all above requirements by system packages manager or compilte them yourself.

    * Compiling toolchain.

        1. GCC >= 6.1 (need support stdc++11)
        2. CMake >= 3.1 (optional, if you using CMake to build this packages)

        On Linux system, we using gcc to compile GCTA, other compile toolchain may works too, but we don't test yet.

        If you want link gcta staticlly, When you use GCC for example, you need all static libs, include stdc++, pthread and gomp. those static lib is part of GCC and GLibc, and you may need compile GCC and GLibc yourself to get those static libs.

2. Compilation
    Here are three different ways to build GCTA, you can chose any one you familiar with to build this package.

    * Using CMake buildsystem generator.

        1. Some dependencies will be installed by CMake FetchContent. Use your system package manager (apt,dnf,pacman,etc.) install the missing ones.
           For example, you may need `sudo pacman -S gsl sqlite3 intel-mkl`. You can also install `ninja` for the `-G Ninja` part in the following CMake command, otherwise just remove `-G Ninja` to fallback to the default (Makefile).

        2. Generate CMake files by:
           `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/build/Release/installed/usr -G Ninja -B build/Release -S .`

        3. Build
           `cmake --build build/Release`

        4. Install
           `cmake --install build/Release`

           Then you can get the executable `build/Release/installed/usr/bin/gcta64`.

        (optional)

        5. You can also package a more portable binary using [linuxdepoly](https://github.com/linuxdeploy/linuxdeploy) and [appimagetool](https://github.com/AppImage/appimagetool). You can install them directly from their `release` page.

            - Generate an AppDir
            `linuxdeploy-x86_64.AppImage --appdir build/Release/installed --executable build/Release/installed/usr/bin/gcta64 --desktop-file build/Release/installed/usr/share/applications/gcta64.desktop --icon-file build/Release/installed/usr/share/icons/hicolor/256x256/apps/gcta64.png`

            - Create an AppImage
            `appimagetool-x86_64.AppImage  build/Release/installed`

            Then you get a portable binary, and you can run it in _almost_ everywhere. Since `AppImage` does not bundle the `glibc` library, you may still run into a incompatible `glibc` problem, if so, you need to build the binary from scratch on that system.

    * Using Makefile (Not implemented yet)

    * Using GNU autoconf/automake system (Not implemented yet)

## macOS

1. Requirments
    * libs
        All libs is as same as libs needed by Linux system. You may need install all libs by yourself or using system's packages manager.

    * compiling toolchain
        1. Clang/LLVM
        2. CMake

        Under macOS, we use Clang/LLVM toolchain to build GCTA, GCC not test yet. Your can using xcode to install Clang/LLVM.

        If you need link GCTA staticlly, you maybe need to compile Clang/LLVM yourself, and LLVM has omp and C standed Lib implement, you may need too compile them too.

2. Complilation

    Procedure is same with Linux system.

## Windows

1. Requirements
    * libs
        All libs is as same as libs needed by Linux system. For intel MKL, Eigen and Spectra, you may need install them yourself. For other libs I recommend you installed them by [vcpkg](https://vcpkg.io/en/index.html). For example to install gsl.
        ```
        vcpkg install gsl:x64-windows-static
        ```
        This command will install gsl static library for 64 bits windows system.

    * Compliling toolchain
        Here, we using VS/Clang toolchain, you can install it by installing Visual Studio. When you install Visual Stedio, you can choose its component, you need choose "Desktop development with c++" workload and "c++ clang compilter for windows" component.

        The Visual Studio(VS for short) is an IDE developed by Microsoft. It contain many components for different languages. MSVC is Microsoft's C/C++ compilter. and which is used internally by VS. And recently, VS have implement Clang/LLVM, and you can use it as VS' C/C++ compiler, or just use it in command line.

        We test MinGW/GCC, which have problem to link MKL. We also test VS/MSVC compilter, which can not compile some source code(plink-ng) because those code using GCC builtin funcitons.

        Clang using most like GCC, the argument is same mostly. you need open VS command prompt, and using it.

2. Compilation

    Until now, due to we have not implement any build system yet, you need to build the GCTA manually.

    First you need compile all source, get object files. and you alse need to compile submods, including `Pgenlib/PgenReader.cpp, plink-ng/2.0/pgenlib_write.cc, plink-ng/2.0/pgenlib_read.cc plink-ng/2.0/pgenlib_misc.cc plink-ng/2.0/pgenlib_ffi_support.cc plink-ng/2.0/plink2_base.cc`.

    Second, link object files and libs, you need link `zlib(Link zlib.lib file. This is different from linux), zstd, gsl, gslblas, sqlite3, mkl_core, mkl_intel_lp64, mkl_intel_thread, libiomp5md/libomp(Link libiomp5md.lib/libomp.lib)`.

    Note that, `libiomp5md` is part of intel compiler, and Clang/LLVM may contain this library too. Of cause you can use `libomp` which is omp implement of Clang/LLVM instead of libiomp5md.

    The dynamic library of windows is different from Unix like system. The `.lib` only contain information for linker, and the code of library is contained by `.dll` file. After link, you may need copy DLL file to the directory which gcta binary located, and when you distribute your binary you need distribute this DLL file with gcta binary too.
