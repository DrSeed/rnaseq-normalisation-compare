# Self-contained demo: a few hugely induced genes make everything else look "down".
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from normalise import METHODS, tpm


def simulate(n_genes=3000, n_per_group=3, n_up=150, fold=10.0, seed=3):
    rng = np.random.default_rng(seed)
    base = rng.lognormal(mean=4.0, sigma=1.6, size=n_genes)
    lengths = rng.lognormal(mean=7.5, sigma=0.6, size=n_genes)
    group = np.array(["ctrl"] * n_per_group + ["treat"] * n_per_group)
    up = np.zeros(n_genes, dtype=bool)
    up[rng.choice(np.argsort(base)[-1000:], n_up, replace=False)] = True
    mu = np.tile(base[:, None], (1, len(group)))
    mu[np.ix_(up, group == "treat")] *= fold
    depth = rng.uniform(0.6, 1.4, len(group))
    mu = mu * depth / mu.sum(axis=0) * 2e7
    disp = 0.05
    counts = rng.negative_binomial(1 / disp, 1 / (1 + mu * disp))
    return counts, lengths, group, up


def log2fc(norm, group):
    a = norm[:, group == "treat"].mean(axis=1) + 0.5
    b = norm[:, group == "ctrl"].mean(axis=1) + 0.5
    return np.log2(a / b)


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    counts, lengths, group, up = simulate()
    expressed = counts.sum(axis=1) >= 20
    lfcs = {name: log2fc(counts / fn(counts), group) for name, fn in METHODS.items()}
    lfcs["TPM"] = log2fc(tpm(counts, lengths), group)

    rows = []
    for name, l in lfcs.items():
        null = l[~up & expressed]
        rows.append({"method": name, "median_lfc_unchanged_genes": round(float(np.median(null)), 3),
                     "pct_unchanged_called_abs_lfc_gt_0.5": round(100 * float(np.mean(np.abs(null) > 0.5)), 1),
                     "median_lfc_induced_genes(true=3.32)": round(float(np.median(l[up])), 2)})
    table = pd.DataFrame(rows)
    table.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    bins = np.linspace(-2, 2, 81)
    for name, l in lfcs.items():
        axes[0].hist(l[~up & expressed], bins=bins, histtype="step", lw=1.6, label=name)
    axes[0].axvline(0, color="k", lw=1)
    axes[0].set_xlabel("log2 fold change of UNCHANGED genes")
    axes[0].set_title("Should be centred on zero")
    axes[0].legend(fontsize=8)

    mean_expr = np.log2(counts.mean(axis=1) + 1)
    for ax, name in zip(axes[1:], ["CPM (library size)", "Median of ratios (DESeq2)"]):
        ax.scatter(mean_expr[~up], lfcs[name][~up], s=3, alpha=0.3, color="grey", label="unchanged")
        ax.scatter(mean_expr[up], lfcs[name][up], s=6, color="#c0504d", label="induced 10x")
        ax.axhline(0, color="k", lw=1)
        ax.set_ylim(-3, 5)
        ax.set_xlabel("mean log2 count")
        ax.set_ylabel("log2 fold change")
        ax.set_title(f"MA plot: {name}")
        ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
