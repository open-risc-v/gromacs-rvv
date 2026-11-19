# 性能与正确性日志 —— `sci-chemistry/gromacs`

每条记录必须包含：**日期 / 变更 / 内核 / 输入规模 / A / B / 比值 / 是否位精确 / 证据路径**。
没有 A/B 双数字的记录不算数。

## 基线（G2，内核关闭）

| 日期 | 构建 | 用例 | 规模 | A(ms) | B(ms) | 备注 |
|---|---|---|---|---|---|---|
|  | riscv64 |  |  |  |  |  |

```
cpu_probe（目标机）：
vlen=
vle_ticks=  vlse_ticks=  vrgather_ticks=  vsetvl_ticks=
```

## 逐轮记录（G3）

| 轮 | 日期 | 内核 | LMUL | 规模 | A(ms) | B(ms) | 收益 | 位精确 | 备注 |
|---|---|---|---|---|---|---|---|---|---|
| 1 |  |  |  |  |  |  |  |  |  |

## 收敛分析（G4）

```
perf top（剩余热点，按占比排序）：
1.
2.
3.
结论：覆盖度 f=____；停止理由：____
```

## 上游构建修复（真机 rv64gcv）

- 失败：`GMXRC.cmakein` / `GMXRC.bash.cmakein` 不存在（`scripts/CMakeLists.txt:39` configure_file）。
- 根因：`src/scripts` 只有 6 个文件（`scripts/GMXRC*` 全缺），`baseline/scripts` 有 10 个；又是复制 src 时被 `.gitignore`/稀疏拷贝吞掉。
- 处置：`rsync -a --ignore-existing baseline/ src/` 补齐到 10 个文件，baseline 未动。
- 结果：`00_framework` 2 文件、`01_gromacs` 2 文件逐字节一致，A/B **ALL PASS**。
