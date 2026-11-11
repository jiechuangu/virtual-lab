"""溶液混合模拟引擎：酸碱平衡 + 沉淀溶解平衡 + 中和热。

模型（25 °C，稀溶液，活度≈浓度）：
- 强电解质完全电离；CH3COOH（Ka=1.75e-5）、NH3·H2O（Kb=1.8e-5）为弱电解质
- AgCl、BaSO4 按 Ksp 精确求解沉淀量；Cu(OH)2 按 Ksp 与 pH 联立求解
- pH 由电荷守恒方程用二分法数值求解
- 温升按中和热 57.3 kJ/mol、溶液比热 4.18 J/(g·K)、密度 1 g/mL 估算
未考虑：配合物（如 [Cu(NH3)4]2+、[Ag(NH3)2]+）、Ag2O 沉淀、离子强度修正、沉淀热。
"""
from dataclasses import dataclass, field
import math

KW = 1.0e-14
KA_HAC = 1.75e-5
KA_NH4 = KW / 1.8e-5
KSP = {"AgCl": 1.8e-10, "BaSO4": 1.1e-10, "Cu(OH)2": 2.2e-20}
DH_NEUT = 57.3e3  # J/mol

# 试剂 -> 每 mol 溶质贡献的物种
REAGENTS = {
    "盐酸 HCl": {"H": 1, "Cl": 1},
    "氢氧化钠 NaOH": {"Na": 1, "OH": 1},
    "醋酸 CH3COOH": {"HAc": 1},
    "氨水 NH3·H2O": {"NH3": 1},
    "氯化钠 NaCl": {"Na": 1, "Cl": 1},
    "硝酸银 AgNO3": {"Ag": 1, "NO3": 1},
    "氯化钡 BaCl2": {"Ba": 1, "Cl": 2},
    "硫酸钠 Na2SO4": {"Na": 2, "SO4": 1},
    "硫酸铜 CuSO4": {"Cu": 1, "SO4": 1},
    "纯水 H2O": {},
}

PRECIP_COLOR = {"AgCl": "#f4f4f4", "BaSO4": "#ffffff", "Cu(OH)2": "#3d7fd6"}
PRECIP_NAME = {"AgCl": "氯化银（白色）", "BaSO4": "硫酸钡（白色）", "Cu(OH)2": "氢氧化铜（蓝色絮状）"}


@dataclass
class Portion:
    reagent: str
    conc: float  # mol/L
    volume_ml: float


@dataclass
class Result:
    volume_ml: float
    pH: float
    temp_c: float
    precipitates: dict  # 名称 -> mol
    free: dict  # 物种 -> mol/L
    reactions: list = field(default_factory=list)


def _pair_precip(a: float, b: float, ksp: float) -> float:
    """1:1 沉淀 AB：返回沉淀浓度 x，使 (a-x)(b-x)=Ksp。"""
    if a * b <= ksp:
        return 0.0
    return ((a + b) - math.sqrt((a - b) ** 2 + 4 * ksp)) / 2


def mix(portions: list[Portion], t0: float = 25.0) -> Result:
    v_l = sum(p.volume_ml for p in portions) / 1000
    if v_l <= 0:
        raise ValueError("总体积必须大于 0")
    tot = {}
    for p in portions:
        for sp, k in REAGENTS[p.reagent].items():
            tot[sp] = tot.get(sp, 0.0) + k * p.conc * p.volume_ml / 1000
    c = {k: v / v_l for k, v in tot.items()}
    g = lambda k: c.get(k, 0.0)

    reactions, prec = [], {}
    # 与 pH 无关的沉淀
    x = _pair_precip(g("Ag"), g("Cl"), KSP["AgCl"])
    ag, cl = g("Ag") - x, g("Cl") - x
    if x > 0:
        prec["AgCl"] = x * v_l
        reactions.append("Ag⁺ + Cl⁻ → AgCl↓")
    y = _pair_precip(g("Ba"), g("SO4"), KSP["BaSO4"])
    ba, so4 = g("Ba") - y, g("SO4") - y
    if y > 0:
        prec["BaSO4"] = y * v_l
        reactions.append("Ba²⁺ + SO₄²⁻ → BaSO₄↓")

    cu_tot, na, no3 = g("Cu"), g("Na"), g("NO3")
    hac, nh3 = g("HAc"), g("NH3")
    # HCl 的 H 与 NaOH 的 OH 不单独出现在电荷守恒中：它们由 Cl⁻/Na⁺ 体现（质子条件）

    def cu_free(h):
        oh = KW / h
        return min(cu_tot, KSP["Cu(OH)2"] / oh**2)

    def f(h):
        pos = h + na + 2 * ba + ag + 2 * cu_free(h) + nh3 * h / (h + KA_NH4)
        neg = KW / h + cl + no3 + 2 * so4 + hac * KA_HAC / (KA_HAC + h)
        return pos - neg

    lo, hi = -16.0, 1.5  # log10[H+]
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(10**mid) > 0:
            hi = mid
        else:
            lo = mid
    h = 10 ** ((lo + hi) / 2)
    pH = -math.log10(h)

    cu = cu_free(h)
    if cu_tot - cu > max(1e-6, 0.01 * cu_tot):  # 忽略水解平衡边缘的痕量值
        prec["Cu(OH)2"] = (cu_tot - cu) * v_l
        reactions.append("Cu²⁺ + 2OH⁻ → Cu(OH)₂↓")

    acid = tot.get("H", 0) + tot.get("HAc", 0)
    base = tot.get("OH", 0) + tot.get("NH3", 0) - 2 * prec.get("Cu(OH)2", 0)
    n_neut = max(0.0, min(acid, base))
    if n_neut > 0:
        reactions.insert(0, "H⁺ + OH⁻ → H₂O（中和反应，放热）")
    dT = n_neut * DH_NEUT / (v_l * 1000 * 4.18)

    free = {
        "H⁺": h, "OH⁻": KW / h, "Na⁺": na, "Cl⁻": cl, "Ag⁺": ag, "Ba²⁺": ba,
        "SO₄²⁻": so4, "Cu²⁺": cu, "NO₃⁻": no3,
        "CH₃COO⁻": hac * KA_HAC / (KA_HAC + h), "CH₃COOH": hac * h / (KA_HAC + h),
        "NH₄⁺": nh3 * h / (h + KA_NH4), "NH₃": nh3 * KA_NH4 / (h + KA_NH4),
    }
    free = {k: v for k, v in free.items() if v > 1e-12}
    return Result(v_l * 1000, pH, t0 + dT, prec, free, reactions)


# 通用指示剂颜色（pH -> RGB）
_UNIV = [(0, (220, 30, 40)), (3, (240, 80, 40)), (4, (250, 150, 40)), (5, (250, 210, 50)),
         (6, (190, 220, 60)), (7, (60, 170, 70)), (8, (40, 150, 140)), (9, (40, 100, 200)),
         (10, (60, 60, 180)), (12, (110, 50, 160)), (14, (90, 30, 120))]


def indicator_rgb(pH: float) -> tuple:
    pH = max(0.0, min(14.0, pH))
    for (p1, c1), (p2, c2) in zip(_UNIV, _UNIV[1:]):
        if pH <= p2:
            t = (pH - p1) / (p2 - p1)
            return tuple(round(a + (b - a) * t) for a, b in zip(c1, c2))
    return _UNIV[-1][1]
