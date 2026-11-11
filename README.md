# Virtual Lab 虚拟实验室

在电脑上做计算化学 / 药物 / 低空经济的模拟实验（in silico）。

> 模拟用于预测与筛选，不能替代真实实验验证。

## 快速开始

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows
pip install -r requirements.txt
python chem/molecules/drug_properties.py
```

## 目录

| 路径 | 内容 |
|---|---|
| `chem/molecules/` | 分子建模与性质预测（RDKit） |
| `chem/docking/` | 药物-靶点对接（计划中） |
| `chem/md/` | 分子动力学（计划中） |
| `aero/` | 低空经济：无人机仿真、航线规划（计划中） |
| `notebooks/` | Jupyter 实验记录本 |
| `results/` | 实验输出数据与图表 |

## 实验日志

| 编号 | 实验 | 脚本 |
|---|---|---|
| 001 | 常见药物分子性质与 Lipinski 类药性分析 | `chem/molecules/drug_properties.py` |
