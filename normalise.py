# RNA-seq count normalisation methods, written out so you can see what each assumes.
import numpy as np


def cpm(counts):
    return counts / counts.sum(axis=0) * 1e6


def tpm(counts, lengths_bp):
    rate = counts / (np.asarray(lengths_bp)[:, None] / 1000.0)
    return rate / rate.sum(axis=0) * 1e6


def _geo_centre(x):
    return x / np.exp(np.mean(np.log(x)))


def size_factors_library(counts):
    return _geo_centre(counts.sum(axis=0).astype(float))


def size_factors_upper_quartile(counts):
    uq = np.array([np.percentile(c[c > 0], 75) for c in counts.T])
    return _geo_centre(uq)


def size_factors_median_ratio(counts):
    # DESeq2: per-sample median ratio to the gene-wise geometric mean.
    with np.errstate(divide="ignore"):
        logc = np.log(counts)
    logmean = logc.mean(axis=1)
    ok = np.isfinite(logmean)
    return np.exp(np.median(logc[ok] - logmean[ok, None], axis=0))


def tmm_factor(obs, ref, trim_m=0.3, trim_a=0.05):
    # edgeR: weighted trimmed mean of log-ratios against a reference sample.
    keep = (obs > 0) & (ref > 0)
    n_o, n_r = obs.sum(), ref.sum()
    o, r = obs[keep], ref[keep]
    m = np.log2((o / n_o) / (r / n_r))
    a = 0.5 * np.log2((o / n_o) * (r / n_r))
    lo_m, hi_m = np.quantile(m, [trim_m, 1 - trim_m])
    lo_a, hi_a = np.quantile(a, [trim_a, 1 - trim_a])
    sel = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
    w = 1.0 / ((n_o - o) / (n_o * o) + (n_r - r) / (n_r * r))
    return 2 ** (np.sum(w[sel] * m[sel]) / np.sum(w[sel]))


def size_factors_tmm(counts):
    props = counts / counts.sum(axis=0)
    uq = np.array([np.percentile(p, 75) for p in props.T])
    ref = int(np.argmin(np.abs(uq - uq.mean())))
    f = _geo_centre(np.array([tmm_factor(counts[:, j], counts[:, ref]) for j in range(counts.shape[1])]))
    return _geo_centre(counts.sum(axis=0) * f)


METHODS = {
    "CPM (library size)": size_factors_library,
    "Upper quartile": size_factors_upper_quartile,
    "TMM (edgeR)": size_factors_tmm,
    "Median of ratios (DESeq2)": size_factors_median_ratio,
}
