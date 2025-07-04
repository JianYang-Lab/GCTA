#!/usr/bin/env sh

# You should cd to the project root
CWD=$(pwd)

# You should set the linuxdeploy and appimagetool path.
echo "Current working directory: ${CWD}"
echo "Path of linuxdeploy ${LINUX_DEPLOY_BIN}"
echo "Path of appimagetool ${APP_IMAGE_TOOL_BIN}"

module load cmake gcc/9.4.0 ninja intelmkl gsl

# Generate cmake
cmake -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=cmake/hpc-toolchain.cmake \
    -DCMAKE_INSTALL_PREFIX=${CWD}/build/Release/installed/usr \
    -G Ninja -B build/Release -S .

# Build
cmake --build build/Release

# Install
cmake --install build/Release

# Packaging
${LINUX_DEPLOY_BIN} --appdir build/Release/installed \
    --executable build/Release/installed/usr/bin/gcta64 \
    --desktop-file build/Release/installed/usr/share/applications/gcta64.desktop \
    --icon-file build/Release/installed/usr/share/icons/hicolor/256x256/apps/gcta64.png

# NOTE: appimagetool need to access network
${APP_IMAGE_TOOL_BIN} build/Release/installed
