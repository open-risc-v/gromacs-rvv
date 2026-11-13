#!/usr/bin/env bash
# 交叉编译基线：容器内把 baseline 编成 riscv64（gmx CLI），qemu-user 跑 editconf，比对哈希。
set -euo pipefail
cd "$(dirname "$0")/.."
IMAGE=${IMAGE:-rvv-cross:2}
MARCH=${MARCH:-rv64gc}
OUT=build/riscv64-b

docker run --rm -v "$PWD":/work -w /work "${IMAGE}" bash -lc "
set -e
rm -rf ${OUT}
mkdir -p ${OUT}
HOST=riscv64-linux-gnu MARCH=${MARCH} \
./cmake/pkg-build.sh /work/baseline /work/${OUT}/pkg > /tmp/pkg.log 2>&1 || { tail -30 /tmp/pkg.log; exit 1; }

mkdir -p ${OUT}/out
qemu-riscv64 -L /usr/riscv64-linux-gnu ${OUT}/pkg/prefix/bin/gmx editconf \
  -f /work/test/data/tiny.gro -o ${OUT}/out/edited.gro -box 5 5 5 -c > /dev/null 2>&1 \
  || { echo 'qemu 运行失败'; exit 1; }
echo \"edited_bytes=\$(wc -c < ${OUT}/out/edited.gro | tr -d ' ')\" > ${OUT}/out/hash.txt
echo \"edited_fnv=\$(sha256sum ${OUT}/out/edited.gro | cut -c1-16)\" >> ${OUT}/out/hash.txt
echo '--- riscv64（${MARCH}）---'; cat ${OUT}/out/hash.txt
echo '--- native 基线 ---'; cat build/ab-data/out/01_gromacs/a/hash.txt
"
