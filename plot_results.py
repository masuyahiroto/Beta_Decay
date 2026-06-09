import matplotlib.pyplot as plt
import numpy as np

DATFILE = "fortran/results.dat"

labels, T_calc, T_exp, logft = [], [], [], []

with open(DATFILE) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        label   = cols[0]
        t_c     = cols[4]
        t_e     = float(cols[5])
        lft     = cols[7]
        if t_c == "NaN" or t_e < 0:
            continue
        labels.append(label)
        T_calc.append(float(t_c))
        T_exp.append(t_e)
        logft.append(float(lft))

T_calc = np.array(T_calc)
T_exp  = np.array(T_exp)
logft  = np.array(logft)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# --- 左: T_calc vs T_exp の比較（棒グラフ） ---
x = np.arange(len(labels))
w = 0.35
axes[0].bar(x - w/2, T_calc, w, label="Calc", color="steelblue")
axes[0].bar(x + w/2, T_exp,  w, label="Exp",  color="tomato")
axes[0].set_yscale("log")
axes[0].set_xticks(x)
axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
axes[0].set_ylabel(r"$T_{1/2}$ [s]")
axes[0].set_title(r"Half-life: Calc vs Exp")
axes[0].legend()
axes[0].grid(axis="y", linestyle="--", alpha=0.5)

# --- 右: T_calc vs T_exp の散布図（対角線が一致） ---
all_t = np.concatenate([T_calc, T_exp])
lim = (all_t.min() * 0.5, all_t.max() * 2.0)
axes[1].plot(lim, lim, "k--", lw=1, label="Calc = Exp")
axes[1].scatter(T_exp, T_calc, color="steelblue", zorder=5)
for i, lb in enumerate(labels):
    axes[1].annotate(lb, (T_exp[i], T_calc[i]),
                     textcoords="offset points", xytext=(5, 3), fontsize=8)
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlim(lim)
axes[1].set_ylim(lim)
axes[1].set_xlabel(r"$T_{1/2}^{\rm exp}$ [s]")
axes[1].set_ylabel(r"$T_{1/2}^{\rm calc}$ [s]")
axes[1].set_title("Calc vs Exp (log scale)")
axes[1].legend()
axes[1].grid(True, linestyle="--", alpha=0.4)

fig.tight_layout()
plt.savefig("results.pdf", bbox_inches="tight")
plt.savefig("results.png", dpi=150, bbox_inches="tight")
print("Saved: results.pdf / results.png")
plt.show()
