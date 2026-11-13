#!/usr/bin/env bash
# 真实包用例：gmx editconf 处理固定 .gro，输出结构文件 + 哈希。
set -euo pipefail
: "${STAGE:=run}" "${BIN:=build/native}" "${DATA:=build/ab-data}" "${OUT:=build/ab-data/out}"
mkdir -p "$DATA" "$OUT"
case "$STAGE" in
  gen) cp "$(dirname "$0")/../data/tiny.gro" "$DATA/tiny.gro" ;;
  run)
    # BIN 是 cmake 构建目录，上游前缀在 <build>/pkg/prefix（不是 BIN 本身）
    "$(dirname "$0")/../pkg_gromacs.sh" "$BIN/pkg/prefix" "$DATA/tiny.gro" "$OUT/edited.gro"
    { echo "edited_bytes=$(wc -c < "$OUT/edited.gro" | tr -d ' ')"; \
      echo "edited_fnv=$(shasum -a 256 "$OUT/edited.gro" 2>/dev/null | cut -c1-16 || sha256sum "$OUT/edited.gro" | cut -c1-16)"; } > "$OUT/hash.txt"
    ;;
  *) echo "unknown STAGE=$STAGE" >&2; exit 2 ;;
esac
