"""
β⁻ 崩壊半減期 計算結果の可視化スクリプト
  入力: fortran/results.dat
  出力: data/results_<timestamp>.dat  (アーカイブ)
        results.pdf / results.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import shutil
import os
from datetime import datetime

# ── パス設定 ────────────────────────────────────────────────
DATFILE  = "fortran/results.dat"
DATADIR  = "data"
OUTPDF   = "results.pdf"
OUTPNG   = "results.png"

os.makedirs(DATADIR, exist_ok=True)

# ── data/ へ自動アーカイブ ─────────────────────────────────
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
archive = os.path.join(DATADIR, f"results_{ts}.dat")
shutil.copy(DATFILE, archive)
print(f"Archived: {archive}")

# ── データ読み込み ─────────────────────────────────────────
labels, Z_list, A_list, Q_list = [], [], [], []
T_calc_all, T_exp_all, logft_all = [], [], []
has_exp = []

with open(DATFILE) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 8:
            continue
        label = cols[0]
        A     = int(cols[1])
        Z     = int(cols[2])
        Q     = float(cols[3])
        t_c   = cols[4]
        t_e   = float(cols[5])
        lft   = cols[7]

        if t_c in ("NaN", "---") or lft == "---":
            continue

        labels.append(label)
        Z_list.append(Z)
        A_list.append(A)
        Q_list.append(Q)
        T_calc_all.append(float(t_c))
        T_exp_all.append(t_e if t_e > 0 else np.nan)
        logft_all.append(float(lft))
        has_exp.append(t_e > 0)

T_calc = np.array(T_calc_all)
T_exp  = np.array(T_exp_all)
logft  = np.array(logft_all)
has_exp = np.array(has_exp)
n = len(labels)
x = np.arange(n)

# ── 比 (実験値がある核種のみ) ──────────────────────────────
ratio = np.where(has_exp, T_calc / T_exp, np.nan)

# ── 色: ratio に応じて緑↔赤 ───────────────────────────────
def ratio_color(r):
    if np.isnan(r):
        return "gray"
    if 0.5 <= r <= 2.0:
        return "#2ecc71"   # 緑: 2倍以内
    if 0.1 <= r <= 10.0:
        return "#f39c12"   # 橙: 10倍以内
    return "#e74c3c"       # 赤: それ以上

colors = [ratio_color(r) for r in ratio]

# ────────────────────────────────────────────────────────────
# 4パネル図
# ────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 11))
fig.suptitle(
    r"$\beta^-$ Decay Half-life: N=126/125/127 Isotopes"
    "\n" r"Phys. Rev. C 109, 064319 (2024)  |  quenching $q=0.74$",
    fontsize=13, fontweight="bold", y=0.98
)

gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.35,
                       left=0.08, right=0.97, top=0.90, bottom=0.10)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

w = 0.35

# ── パネル1: T_calc vs T_exp 棒グラフ ─────────────────────
bars_c = ax1.bar(x - w/2, T_calc, w, label="Calc", color="steelblue", alpha=0.85)
bars_e = ax1.bar(x + w/2, T_exp,  w, label="Exp",  color="tomato",    alpha=0.85)
ax1.set_yscale("log")
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
ax1.set_ylabel(r"$T_{1/2}$ [s]", fontsize=10)
ax1.set_title("(1) Half-life: Calc vs Exp", fontsize=10, fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(axis="y", linestyle="--", alpha=0.4)
ax1.set_ylim(bottom=0.1)

# ── パネル2: Calc/Exp 比 ──────────────────────────────────
mask = has_exp
ax2.bar(x[mask], ratio[mask], color=[colors[i] for i in range(n) if has_exp[i]],
        alpha=0.85, edgecolor="white")
ax2.axhline(1.0, color="black",   lw=1.5, ls="-",  label="Calc = Exp")
ax2.axhline(2.0, color="#2ecc71", lw=1.0, ls="--", label="×2 limit")
ax2.axhline(0.5, color="#2ecc71", lw=1.0, ls="--")
ax2.axhline(10., color="#f39c12", lw=1.0, ls=":",  label="×10 limit")
ax2.axhline(0.1, color="#f39c12", lw=1.0, ls=":")
ax2.set_yscale("log")
ax2.set_xticks(x[mask])
ax2.set_xticklabels([labels[i] for i in range(n) if has_exp[i]],
                    rotation=35, ha="right", fontsize=8)
ax2.set_ylabel(r"$T_{\rm calc} / T_{\rm exp}$", fontsize=10)
ax2.set_title("(2) Ratio Calc/Exp  (green ≤ ×2)", fontsize=10, fontweight="bold")
ax2.legend(fontsize=8, loc="upper right")
ax2.grid(axis="y", linestyle="--", alpha=0.4)

# ratio の数値をバーの上に表示
for i in range(n):
    if has_exp[i] and not np.isnan(ratio[i]):
        yi = ratio[i]
        va = "bottom" if yi >= 1 else "top"
        offset = 1.15 if yi >= 1 else 0.85
        ax2.text(x[i], yi * offset, f"{yi:.1f}",
                 ha="center", va=va, fontsize=7.5, fontweight="bold")

# ── パネル3: 散布図 T_exp vs T_calc ──────────────────────
mask_sc = has_exp & ~np.isnan(T_exp)
T_c_sc = T_calc[mask_sc]
T_e_sc = T_exp[mask_sc]
labels_sc = [labels[i] for i in range(n) if has_exp[i] and not np.isnan(T_exp[i])]
colors_sc  = [colors[i] for i in range(n) if has_exp[i] and not np.isnan(T_exp[i])]

all_t = np.concatenate([T_c_sc, T_e_sc])
lim = (all_t.min() * 0.3, all_t.max() * 3.0)

ax3.fill_between(lim, [l*0.5 for l in lim], [l*2.0 for l in lim],
                 alpha=0.12, color="#2ecc71", label="×2 band")
ax3.fill_between(lim, [l*0.1 for l in lim], [l*10. for l in lim],
                 alpha=0.08, color="#f39c12", label="×10 band")
ax3.plot(lim, lim, "k-",  lw=1.5, label="Calc = Exp")
ax3.plot(lim, [l*2   for l in lim], "g--", lw=0.8, alpha=0.6)
ax3.plot(lim, [l*0.5 for l in lim], "g--", lw=0.8, alpha=0.6)
ax3.scatter(T_e_sc, T_c_sc, c=colors_sc, s=70, zorder=5, edgecolors="gray", lw=0.5)
for i, lb in enumerate(labels_sc):
    ax3.annotate(lb, (T_e_sc[i], T_c_sc[i]),
                 textcoords="offset points", xytext=(6, 3),
                 fontsize=7.5, color="dimgray")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(lim)
ax3.set_ylim(lim)
ax3.set_xlabel(r"$T_{1/2}^{\rm exp}$ [s]", fontsize=10)
ax3.set_ylabel(r"$T_{1/2}^{\rm calc}$ [s]", fontsize=10)
ax3.set_title("(3) Scatter: Calc vs Exp (log–log)", fontsize=10, fontweight="bold")
ax3.legend(fontsize=8, loc="upper left")
ax3.grid(True, linestyle="--", alpha=0.35)

# ── パネル4: log ft ────────────────────────────────────────
bar_colors4 = ["#3498db" if not np.isnan(lft) else "gray" for lft in logft]
ax4.bar(x, logft, color=bar_colors4, alpha=0.85, edgecolor="white")
ax4.axhspan(4.0, 6.0, alpha=0.08, color="green",  label="Allowed GT (4–6)")
ax4.axhspan(6.0, 9.0, alpha=0.08, color="orange", label="First forbidden (6–9)")
ax4.axhline(5.0, color="green",  lw=1.0, ls="--", alpha=0.7)
ax4.axhline(7.0, color="orange", lw=1.0, ls="--", alpha=0.7)
ax4.set_xticks(x)
ax4.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
ax4.set_ylabel(r"$\log_{10} ft$", fontsize=10)
ax4.set_title("(4) log ft  (placeholder B_GT)", fontsize=10, fontweight="bold")
ax4.legend(fontsize=8, loc="upper right")
ax4.grid(axis="y", linestyle="--", alpha=0.4)
ax4.set_ylim(0, max(logft) * 1.15)

for i, lft_v in enumerate(logft):
    ax4.text(x[i], lft_v + 0.05, f"{lft_v:.2f}",
             ha="center", va="bottom", fontsize=7.5)

# ── 凡例パッチ (ratio カラー説明) ─────────────────────────
legend_patches = [
    mpatches.Patch(color="#2ecc71", label=r"ratio ≤ ×2 (good)"),
    mpatches.Patch(color="#f39c12", label=r"ratio ≤ ×10"),
    mpatches.Patch(color="#e74c3c", label=r"ratio > ×10"),
    mpatches.Patch(color="gray",    label="no exp. data"),
]
fig.legend(handles=legend_patches, loc="lower center", ncol=4,
           fontsize=9, framealpha=0.8,
           bbox_to_anchor=(0.5, 0.01))

# ── 保存 ──────────────────────────────────────────────────
plt.savefig(OUTPDF, bbox_inches="tight")
plt.savefig(OUTPNG, dpi=150, bbox_inches="tight")
print(f"Saved: {OUTPDF}")
print(f"Saved: {OUTPNG}")
