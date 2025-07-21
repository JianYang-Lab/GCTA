#!/usr/bin/env sh

set -e
set -o pipefail

# You should cd to the project root
CWD=$(pwd)

# You should set the linuxdeploy and appimagetool path.
echo "Current working directory: ${CWD}"
echo "Path of linuxdeploy ${LINUX_DEPLOY_BIN}"
echo "Path of appimagetool ${APP_IMAGE_TOOL_BIN}"
echo "Path of appimagetool ${APP_IMAGE_RUNTIME_FILE}"

module load gcc/8.4.0 cmake intelmkl gsl

# Generate cmake
cmake -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=cmake/hpc-toolchain.cmake \
    -DCMAKE_INSTALL_PREFIX=${CWD}/build/Release/installed/usr \
    -B build/Release -S .

# Build
cmake --build build/Release

# Install
cmake --install build/Release

cp /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/intel64/libmkl_avx512.so.2 build/Release/installed/usr/lib/
strip build/Release/installed/usr/lib/libmkl_avx512.so.2

# Packaging
${LINUX_DEPLOY_BIN} --appdir build/Release/installed \
    --executable build/Release/installed/usr/bin/gcta64 \
    --desktop-file build/Release/installed/usr/share/applications/gcta64.desktop \
    --icon-file build/Release/installed/usr/share/icons/hicolor/256x256/apps/gcta64.png

# NOTE: appimagetool need to access network
${APP_IMAGE_TOOL_BIN} build/Release/installed --runtime-file ${APP_IMAGE_RUNTIME_FILE}
