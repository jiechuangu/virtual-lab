"""实验 001：常见药物的分子性质与类药性（Lipinski 五规则）分析。

运行：python chem/molecules/drug_properties.py
输出：results/001_drug_properties.csv、results/001_drugs.png、results/001_logp_vs_mw.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Lipinski, QED, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results"

# 名称 -> SMILES（分子的文本表示法）
DRUGS = {
    "Aspirin 阿司匹林": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Ibuprofen 布洛芬": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Paracetamol 对乙酰氨基酚": "CC(=O)NC1=CC=C(O)C=C1",
    "Caffeine 咖啡因": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Metformin 二甲双胍": "CN(C)C(=N)NC(=N)N",
    "Penicillin G 青霉素G": "CC1(C)SC2C(NC(=O)CC3=CC=CC=C3)C(=O)N2C1C(=O)O",
    "Atorvastatin 阿托伐他汀": "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
}


def analyze(name: str, smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析 {name}: {smiles}")
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "name": name,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "MW": round(mw, 1),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "QED": round(QED.qed(mol), 3),
        "Lipinski_violations": violations,
        "drug_like": violations <= 1,
    }


def main() -> None:
    OUT.mkdir(exist_ok=True)
    df = pd.DataFrame([analyze(n, s) for n, s in DRUGS.items()])
    df.to_csv(OUT / "001_drug_properties.csv", index=False, encoding="utf-8-sig")
    print(df.to_string(index=False))

    mols = [Chem.MolFromSmiles(s) for s in DRUGS.values()]
    legends = [n.split(" ")[0] for n in DRUGS]
    Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300), legends=legends).save(
        OUT / "001_drugs.png"
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["MW"], df["LogP"], s=df["QED"] * 400, alpha=0.6)
    for _, r in df.iterrows():
        ax.annotate(r["name"].split(" ")[0], (r["MW"], r["LogP"]), fontsize=8)
    ax.axvline(500, ls="--", c="gray")
    ax.axhline(5, ls="--", c="gray")
    ax.set_xlabel("Molecular weight (Da)")
    ax.set_ylabel("LogP")
    ax.set_title("Lipinski space (bubble size = QED)")
    fig.tight_layout()
    fig.savefig(OUT / "001_logp_vs_mw.png", dpi=150)
    print(f"\n结果已保存到 {OUT}")


if __name__ == "__main__":
    main()
