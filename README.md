# 论文代码与结果

```text
EHD_improve_1115/
├── README.md
├── algorithm/
│   ├── code/          共享模型与算法比较代码
│   ├── results/       算法比较数据及汇总
│   └── figure/        .fig、.png
├── sensitive/
│   ├── code/          参数分析代码
│   ├── results/       参数敏感性数据
│   └── figure/        .fig、.png
└── MG_epsilon/
    ├── voltage/
    │   ├── code/      固定/动态介电及电场力分解代码
    │   ├── results/   不同电压的数据
    │   └── figure/    .fig、.png
    └── rise_time/
        ├── code/      固定/动态介电比较代码
        ├── results/   不同升压时间的数据
        └── figure/    .fig、.png
```

代码为 `.m`，数据为 `.csv`、`.mat`，图件为 `.fig`、`.png`；必要的 `.json` 保存运行状态和绘图配置。

各分类共用 `algorithm/code/model/` 中的物理模型和求解器。完整输入表为 `algorithm/code/case_catalog.csv`；动态两类比较共用的 50 kV、1 ns 数据位于电压比较目录，升压时间比较通过 `results/shared_baseline.csv` 引用。

使用 MATLAB R2025b：进入 `algorithm/code/`，运行 `plot_all` 重绘已有结果，运行 `run_all` 重新计算全部 47 个算例（会覆盖对应输出）。
