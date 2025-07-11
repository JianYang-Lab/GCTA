set(LLVM_ROOT "C:/Program Files/LLVM" CACHE STRING "LLVM install directory")
set(CMAKE_C_COMPILER "${LLVM_ROOT}/bin/clang.exe")
set(CMAKE_CXX_COMPILER "${LLVM_ROOT}/bin/clang++.exe")
# plink compiler version assertion
add_compile_options(-fgnuc-version=40900)

