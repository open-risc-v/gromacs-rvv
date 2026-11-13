#!/usr/bin/env bash
# 构建上游 gromacs（cmake）→ <out>/prefix。
# 基线口径：GMX_FFT_LIBRARY=fftpack（内置，无外部 FFT 依赖）、GMX_MPI=OFF、
# GMX_SIMD=None（aarch64 的 NEON 与 riscv 的标量结果不保证逐位相同，基线先走同一份标量代码）。
set -euo pipefail
SRC=$1
OUT=$2
rm -rf "${OUT}"
mkdir -p "${OUT}/prefix" "${OUT}/src"
cp -R "${SRC}" "${OUT}/src/pkg"
cd "${OUT}/src/pkg"
EXTRA=""
if [ -n "${HOST:-}" ]; then
  EXTRA="-DCMAKE_TOOLCHAIN_FILE=/work/cmake/toolchain-riscv64-linux.cmake -DCMAKE_C_FLAGS=-march=${MARCH:-rv64gc} -DCMAKE_CXX_FLAGS=-march=${MARCH:-rv64gc}"
fi
cmake -S . -B "${OUT}/bld" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${OUT}/prefix" \
  -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=OFF -DGMX_THREAD_MPI=OFF -DGMX_GPU=OFF \
  -DGMX_SIMD=None -DGMX_DOUBLE=OFF -DGMX_BUILD_UNITTESTS=OFF \
  -DBUILD_TESTING=OFF -DGMX_INSTALL_LEGACY_API=OFF ${EXTRA} >/dev/null
ninja -C "${OUT}/bld" >/dev/null
ninja -C "${OUT}/bld" install >/dev/null
echo "== 上游 gromacs 构建完成：${OUT}/prefix"
