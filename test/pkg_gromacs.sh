#!/usr/bin/env bash
# gromacs 探针：grompp + 极短 mdrun（单分子水盒，20 步），比末帧坐标与轨迹字节。
# 选 mdrun 而不是 editconf：editconf 的主路径只有坐标平移与包围盒统计，
# 没有可位精确向量化的浮点热点；SETTLE 约束则在每个积分步都跑，逐分子独立。
# 该构建 GMX_THREAD_MPI=OFF，mdrun 不能传 -nt（会致命报错）。
set -euo pipefail
PREFIX=$1
DATADIR=$2
OUTDIR=$3
export GMXLIB="${PREFIX}/share/gromacs/top"
mkdir -p "${OUTDIR}"
cp "${DATADIR}/tiny-water.gro" "${OUTDIR}/w.gro"
"${PREFIX}/bin/gmx" grompp -f "${DATADIR}/grompp.mdp" -c "${OUTDIR}/w.gro" \
  -p "${DATADIR}/topol.top" -o "${OUTDIR}/topol.tpr" > /dev/null 2>&1
"${PREFIX}/bin/gmx" mdrun -s "${OUTDIR}/topol.tpr" -o "${OUTDIR}/traj.trr" \
  -c "${OUTDIR}/confout.gro" > /dev/null 2>&1
