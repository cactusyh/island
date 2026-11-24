# J101 · 固定 75% 环氧转化率的盒尺寸与采样对照

使用 J97 同一 PCFF 库、450 K/1 atm 固化协议，设定 75% 的反应目标。在三种盒尺寸（4+2、8+4、16+8 个 DGEBA+PACM）及四个装填种子（41、173、2026、719）上分别执行。提前停滞的样品计入总数，只在实际达到相同转化率时进入性质比较。对每个入选样品的固化前、后结构分别进行 300 K/1 atm NPT。两种观察窗口取同一条轨迹前 6 ps 和前 20 ps 各自的后半段，由每个种子的配对平均体积计算收缩率。95% 区间按独立装填种子估计，同时保留平均值、种子间标准差和窗口间配对变化。

```bash
python scripts/benchmark_cure_size.py work/size --library work/library/library.json --workers 3
python -m island.epoxy.convergence work/size/cases.json
python scripts/benchmark_cure_rates.py work/size/d8-p4-s173 work/rates8 --workers 2
python scripts/benchmark_cure_rates.py work/size/d16-p8-s173 work/rates16 --workers 2
```

冷却用 550–250 K、4 和 8 K/ps；拉伸在 300 K、0.001 和 0.002 ps⁻¹ 下各沿三个方向。沿用 J97 的折点及方向模量质量判据。此处树脂为纯 DGEBA n=0 / PACM，文献中的商业树脂密度约 1.2 g/cm³，固化度和热历史有差异。[实验比较与拟合定义](phase_4j95_j97.md)
