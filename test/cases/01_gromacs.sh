#!/usr/bin/env bash
# 真实包用例：grompp + 极短 mdrun，输出末帧坐标与轨迹的字节哈希。
set -euo pipefail
: "${STAGE:=run}" "${SIDE:=A}" "${BIN:=build/native}" "${DATA:=build/ab-data}" "${OUT:=build/ab-data/out}"
mkdir -p "$DATA" "$OUT"
case "$STAGE" in
  gen) cp "$(dirname "$0")/../data/tiny-water.gro" "$DATA/tiny-water.gro"
       cp "$(dirname "$0")/../data/grompp.mdp" "$DATA/grompp.mdp"
       cp "$(dirname "$0")/../data/topol.top" "$DATA/topol.top" ;;
  run)
    "$(dirname "$0")/../pkg_gromacs.sh" "$BIN/pkg/prefix" "$DATA" "$OUT"
    { echo "confout_bytes=$(wc -c < "$OUT/confout.gro" | tr -d ' ')";
      echo "confout_fnv=$( (shasum -a 256 "$OUT/confout.gro" 2>/dev/null || sha256sum "$OUT/confout.gro") | cut -c1-16)";
      echo "traj_bytes=$(wc -c < "$OUT/traj.trr" | tr -d ' ')";
      echo "traj_fnv=$( (shasum -a 256 "$OUT/traj.trr" 2>/dev/null || sha256sum "$OUT/traj.trr") | cut -c1-16)"; } > "$OUT/hash.txt"
    ;;
  *) echo "unknown STAGE=$STAGE" >&2; exit 2 ;;
esac
