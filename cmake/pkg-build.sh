#!/usr/bin/env bash
# 构建上游 gromacs（cmake）→ <out>/prefix。
# 基线口径：GMX_FFT_LIBRARY=fftpack、GMX_MPI=OFF、GMX_SIMD=None（保证 A/B 走同一份标量代码）。
# A/B 差异只由位置参数驱动：$3 = -march 值，$4 = 1 表示打开 RVV 内核。
# 内核必须随上游一起编，所以宏与 -march 通过 CMAKE_C_FLAGS/CXX_FLAGS 传进上游。
set -euo pipefail
SRC=$1
OUT=$2
MARCH=${3:-rv64gc}
KERN=${4:-0}
rm -rf "${OUT}"
mkdir -p "${OUT}/prefix" "${OUT}/src"
OUT=$(realpath "${OUT}")
cp -R "${SRC}" "${OUT}/src/pkg"
cd "${OUT}/src/pkg"
CMF="-march=${MARCH}"
if [ "${KERN}" = "1" ]; then
  CMF="${CMF} -DRVV_KERNELS"
fi
# 注意：-D 选项必须拆成两个 argv，且每个整体加引号；
# 若把两个 -D 拼成一个字符串，cmake 会把第二个当成前一个的值。
EXTRA=(-DCMAKE_C_FLAGS="${CMF}" -DCMAKE_CXX_FLAGS="${CMF}")
cmake -S . -B "${OUT}/bld" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${OUT}/prefix" \
  -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=OFF -DGMX_THREAD_MPI=OFF -DGMX_GPU=OFF \
  -DGMX_SIMD=None -DGMX_DOUBLE=OFF -DGMX_BUILD_UNITTESTS=OFF \
  -DBUILD_TESTING=OFF -DGMX_INSTALL_LEGACY_API=OFF "${EXTRA[@]}" >/dev/null
ninja -C "${OUT}/bld" >/dev/null
ninja -C "${OUT}/bld" install >/dev/null
echo "== 上游 gromacs 构建完成：${OUT}/prefix"
