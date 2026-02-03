#!/bin/bash

set -e

# 项目配置
APP_NAME="gcta64"
BUILD_DIR="build/Release"
INSTALL_DIR="$BUILD_DIR/installed"
PACKAGE_DIR="${APP_NAME}-package"
EXECUTABLE_PATH="$INSTALL_DIR/usr/bin/$APP_NAME"

# 构建和安装
echo "[1/5] CMake 构建..."
cmake -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/usr" \
  -B "$BUILD_DIR" -S .
cmake --build "$BUILD_DIR" --target install

# 准备打包目录
echo "[2/5] 准备打包目录..."
rm -rf "$PACKAGE_DIR"
mkdir -p "$PACKAGE_DIR/bin"

cp "$EXECUTABLE_PATH" "$PACKAGE_DIR/bin/$APP_NAME"

#dylib 路径
echo "[3/5] 修复 dylib 路径..."

dylibbundler -od -b \
  -x "$PACKAGE_DIR/bin/$APP_NAME" \
  -d "$PACKAGE_DIR/lib" \
  -p @executable_path/../lib \
  -s "$INSTALL_DIR/usr/lib"

#生成 run.sh 启动脚本
echo "[4/5] 创建运行脚本..."
cat > "$PACKAGE_DIR/run.sh" <<EOF
#!/bin/bash
DIR=\$(cd "\$(dirname "\$0")"; pwd)
export DYLD_LIBRARY_PATH="\$DIR/lib"
"\$DIR/bin/$APP_NAME" "\$@"
EOF

chmod +x "$PACKAGE_DIR/run.sh"
