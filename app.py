"""虚拟实验室 UI。运行：streamlit run app.py"""
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "chem" / "solution_lab"))
from engine import PRECIP_COLOR, PRECIP_NAME, REAGENTS, Portion, indicator_rgb, mix  # noqa: E402

st.set_page_config(page_title="Virtual Lab 虚拟实验室", page_icon="🧪", layout="wide")
st.title("🧪 溶液混合实验台")
st.caption("基于酸碱平衡、沉淀溶解平衡与中和热的定量模拟（25 °C，稀溶液近似）")

NAMES = list(REAGENTS)


def beaker_svg(res, indicator: bool) -> str:
    # 溶液颜色
    if indicator:
        r, g, b = indicator_rgb(res.pH)
        alpha = 0.75
    else:
        cu = res.free.get("Cu²⁺", 0)
        r, g, b = (80, 160, 230) if cu > 0 else (200, 225, 245)
        alpha = min(0.75, 0.15 + cu * 4) if cu > 0 else 0.25
    level = min(1.0, res.volume_ml / 250)
    top = 40 + 200 * (1 - level)
    liquid_h = 240 - top
    parts = [
        f'<rect x="42" y="{top:.0f}" width="156" height="{liquid_h:.0f}" fill="rgba({r},{g},{b},{alpha})"/>'
    ]
    # 沉淀：浑浊颗粒 + 底部沉积层
    total_prec = sum(res.precipitates.values())
    y_bottom = 240
    for name, mol in res.precipitates.items():
        frac = mol / total_prec
        h = max(4, min(40, mol * 1000 * 6)) * frac
        col = PRECIP_COLOR[name]
        y_bottom -= h
        parts.append(f'<rect x="42" y="{y_bottom:.0f}" width="156" height="{h:.0f}" fill="{col}" opacity="0.95"/>')
        for i in range(min(60, int(mol * 1000 * 15) + 8)):
            x = 50 + (i * 37) % 140
            y = top + 8 + (i * 53) % max(10, int(liquid_h - 20))
            parts.append(f'<circle cx="{x}" cy="{y:.0f}" r="2" fill="{col}" opacity="0.8"/>')
    ticks = "".join(
        f'<line x1="198" y1="{40 + 200 * (1 - v / 250):.0f}" x2="188" y2="{40 + 200 * (1 - v / 250):.0f}" stroke="#666"/>'
        f'<text x="202" y="{44 + 200 * (1 - v / 250):.0f}" font-size="10" fill="#666">{v}</text>'
        for v in (50, 100, 150, 200, 250)
    )
    return f"""
<svg viewBox="0 0 240 260" width="300" xmlns="http://www.w3.org/2000/svg">
  {''.join(parts)}
  <path d="M40 30 L40 242 L200 242 L200 30" fill="none" stroke="#555" stroke-width="3"/>
  {ticks}
  <text x="120" y="20" text-anchor="middle" font-size="12" fill="#444">烧杯 250 mL</text>
</svg>"""


tab_mix, tab_titr = st.tabs(["🧫 混合实验", "📈 滴定曲线"])

with tab_mix:
    left, right = st.columns([1, 1])
    with left:
        st.subheader("加入试剂")
        n = st.number_input("试剂份数", 1, 5, 2)
        portions = []
        defaults = ["盐酸 HCl", "氢氧化钠 NaOH", "硫酸铜 CuSO4", "硝酸银 AgNO3", "纯水 H2O"]
        for i in range(int(n)):
            c1, c2, c3 = st.columns([2, 1, 1])
            name = c1.selectbox(f"试剂 {i + 1}", NAMES, NAMES.index(defaults[i]), key=f"r{i}")
            conc = 0.0 if name == "纯水 H2O" else c2.number_input(
                "浓度 mol/L", 0.0, 2.0, 0.1, 0.01, key=f"c{i}")
            vol = c3.slider("体积 mL", 0, 100, 50, key=f"v{i}")
            portions.append(Portion(name, conc, vol))
        indicator = st.checkbox("滴加通用指示剂（颜色随 pH 变化）", True)

    with right:
        st.subheader("实验结果")
        if sum(p.volume_ml for p in portions) == 0:
            st.info("请加入试剂")
        else:
            res = mix(portions)
            b1, b2 = st.columns([1, 1])
            b1.markdown(beaker_svg(res, indicator), unsafe_allow_html=True)
            b2.metric("pH", f"{res.pH:.2f}")
            b2.metric("温度", f"{res.temp_c:.2f} °C", f"{res.temp_c - 25:+.2f} °C")
            b2.metric("总体积", f"{res.volume_ml:.0f} mL")
            if res.reactions:
                st.markdown("**发生的反应**")
                for rx in res.reactions:
                    st.markdown(f"- {rx}")
            else:
                st.markdown("**未发生明显反应**")
            for name, mol in res.precipitates.items():
                st.success(f"生成沉淀：{PRECIP_NAME[name]} {mol * 1000:.3f} mmol")
            with st.expander("溶液中各离子浓度"):
                st.dataframe(
                    pd.DataFrame({"物种": list(res.free), "浓度 mol/L": [f"{v:.3e}" for v in res.free.values()]}),
                    hide_index=True,
                )

with tab_titr:
    st.subheader("滴定：逐渐加入滴定剂，观察 pH 变化")
    c1, c2 = st.columns(2)
    with c1:
        a = st.selectbox("烧杯中的溶液", NAMES, NAMES.index("醋酸 CH3COOH"), key="ta")
        ac = st.number_input("浓度 mol/L", 0.001, 2.0, 0.1, key="tac")
        av = st.number_input("体积 mL", 1.0, 100.0, 25.0, key="tav")
    with c2:
        b = st.selectbox("滴定剂（滴定管中）", NAMES, NAMES.index("氢氧化钠 NaOH"), key="tb")
        bc = st.number_input("浓度 mol/L", 0.001, 2.0, 0.1, key="tbc")
        bmax = st.number_input("最多滴加 mL", 1.0, 200.0, 50.0, key="tbm")
    vols = [bmax * i / 200 for i in range(201)]
    ph = [mix([Portion(a, ac, av), Portion(b, bc, v)]).pH for v in vols]
    st.line_chart(pd.DataFrame({"滴加体积 mL": vols, "pH": ph}), x="滴加体积 mL", y="pH")
    cur = st.slider("查看某一点：已滴加 mL", 0.0, float(bmax), float(bmax) / 2)
    r = mix([Portion(a, ac, av), Portion(b, bc, cur)])
    st.markdown(f"滴加 **{cur:.1f} mL** 时：pH = **{r.pH:.2f}**，温度 {r.temp_c:.2f} °C")

with st.sidebar:
    st.header("关于")
    st.markdown(
        "- 模型：电荷守恒 + Ka/Kb + Ksp 数值求解\n"
        "- 中和热 57.3 kJ/mol\n"
        "- **未考虑**：配合物（如铜氨深蓝色）、Ag₂O 沉淀、离子强度修正\n\n"
        "模拟结果用于学习与预测，不能替代真实实验。"
    )
