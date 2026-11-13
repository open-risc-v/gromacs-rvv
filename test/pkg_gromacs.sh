#!/usr/bin/env bash
# gromacs 探针：用 gmx editconf 把固定 .gro 加盒子，比输出结构文件。
# 选 editconf 而不是 mdrun/grompp：它走真实的上游代码路径（拓扑/坐标解析、盒子处理），
# 但不需要力场、不需要 .tpr，基线阶段最省事且完全确定。
set -euo pipefail
PREFIX=$1
IN=$2
OUT=$3
GMXLIB="${PREFIX}/share/gromacs/top" "${PREFIX}/bin/gmx" editconf \
  -f "${IN}" -o "${OUT}" -box 5 5 5 -c > /dev/null 2>&1
