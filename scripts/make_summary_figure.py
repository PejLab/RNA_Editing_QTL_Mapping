#!/usr/bin/env python3
# Summary figure of cis-edQTL results across GTEx tissues, in the style of Li et al. 2022 Fig. 1.
# Inputs: results/figures/data/ (from scripts/prepare_figure_data.py), plus for panel a the editing matrix and genotypes
# of one tissue and the site-to-gene annotation.
#   python3 scripts/make_summary_figure.py

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patheffects
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from pandas_plink import read_plink1_bin

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = f"{ROOT}/results/figures/data"
GENE_ANNOT = f"{ROOT}/src/edit_site_gene_annotation.tsv.gz"   # from scripts/annotate_edit_sites.R (optional)
OUT = f"{ROOT}/results/figures/edqtl_summary.png"
EXAMPLE_TISSUE = "Brain_Cerebellum"

BLUE, ORANGE = "#2a78d6", "#eb6834"            # validated categorical slots 1-2: other tissues, brain
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
})


def panel_label(ax, s, x=-0.2):
    ax.text(x, 1.06, s, transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK, va="bottom")


def pretty(t):
    return t.replace("Brain_", "").replace("_basal_ganglia", "").replace("_", " ").replace("Cells ", "")


summary = pd.read_csv(f"{DATA}/tissue_summary.tsv", sep="\t")
lead = pd.read_csv(f"{DATA}/lead_edqtls.tsv.gz", sep="\t")
cross = pd.read_csv(f"{DATA}/cross_tissue_effects.tsv.gz", sep="\t")
summary["brain"] = summary.tissue.str.startswith("Brain_")
lead["brain"] = lead.tissue.str.startswith("Brain_")
col = lambda b: ORANGE if b else BLUE
print(f"{len(summary)} tissues, {summary.samples.sum():,} samples, {summary.edsites.sum():,} edSites")

fig = plt.figure(figsize=(7.2, 5.3))   # Nature double column: 183 mm wide
gs = fig.add_gridspec(2, 3, left=0.08, right=0.985, top=0.93, bottom=0.08, wspace=0.6, hspace=0.6)
ax = {k: fig.add_subplot(gs[i, j]) for k, (i, j) in zip("abcdef", [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])}
group_legend = [Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="Other tissues"),
                Line2D([], [], marker="o", ls="", color=ORANGE, ms=4, label="Brain")]

# ---------------- a: example edQTL (editing level by genotype) ----------------
t = EXAMPLE_TISSUE
c = lead[(lead.tissue == t) & lead.af.between(0.2, 0.8)].sort_values("pval_beta")
mat = pd.read_csv(f"{ROOT}/intermediate/edit_matrix/{t}/{t}.edMat.20cov.60samps.noXYM.txt", sep=" ", index_col=0)
key = lambda pid: f"{pid.rsplit('_', 1)[0]}:{int(pid.rsplit('_', 1)[1]) - 1}:{pid.rsplit('_', 1)[1]}"
row = next(r for _, r in c.iterrows() if key(r.phenotype_id) in mat.index)
site, var = row.phenotype_id, row.variant_id
G = read_plink1_bin(*(f"{ROOT}/data/genotype/{t}/{t}.geno.clean.maf0.05.{e}" for e in ["bed", "bim", "fam"]), verbose=False)
sel = G.snp.values == var
g = pd.Series(G.sel(variant=G.variant[sel]).values.ravel(), index=G.sample.values)
counts = mat.loc[key(site)].str.split("/", expand=True).astype(int)
lvl = (counts[0] / counts[1])[counts[1] >= 20]                    # samples with >= 20 reads
df = pd.DataFrame({"lvl": lvl, "g": g.reindex(lvl.index)}).dropna()
ref, alt = var.split("_")[2:4]
n_alt = (2 - df.g) if G.a0.values[sel][0] == ref else df.g          # plink counts a0 alleles
df["gt"] = n_alt.map({0: f"{ref}{ref}", 1: f"{ref}{alt}", 2: f"{alt}{alt}"})
order = [f"{ref}{ref}", f"{ref}{alt}", f"{alt}{alt}"]
ch, x = site.rsplit("_", 1)
gene = None
if os.path.exists(GENE_ANNOT):                                    # keyed by chrN_<1-based position>, like the site IDs
    annot = pd.read_csv(GENE_ANNOT, sep="\t", usecols=["chr_pos_edsite", "gene_name", "sense"])
    hit = annot[annot.chr_pos_edsite == site].sort_values("sense", ascending=False).gene_name.dropna()   # sense gene first
    gene = hit.iloc[0] if len(hit) else None
a = ax["a"]
rng = np.random.default_rng(1)
for i, gt in enumerate(order):
    v = df.loc[df["gt"] == gt, "lvl"].values
    if len(v):
        a.boxplot(v, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor="#efeeea", edgecolor=INK2, linewidth=0.7),
                  medianprops=dict(color=INK, linewidth=1.1), whiskerprops=dict(color=INK2, linewidth=0.7),
                  capprops=dict(color=INK2, linewidth=0.7))
        a.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=6, color=ORANGE, alpha=0.75, linewidths=0, zorder=3)
ns = [int((df["gt"] == gt).sum()) for gt in order]
a.set_xticks(range(3)); a.set_xticklabels([f"{gt}\n(n = {n})" for gt, n in zip(order, ns)])
vch, vpos = var.split("_")[:2]
a.set_title((f"{gene}, " if gene else "") + f"edited site {ch}:{int(x):,}", color=INK)
a.set_xlabel(f"Genotype at {vch}:{int(vpos):,} ({ref}>{alt})")
a.set_ylabel("Editing level (edited/total reads)")
med = [df.loc[df["gt"] == gt, "lvl"].median() for gt in order if (df["gt"] == gt).any()]
left = med[0] < med[-1]                                             # put the note in the emptier top corner
a.text(0.04 if left else 0.96, 0.97, f"{pretty(t)}\nq = {row.qval:.1e}", transform=a.transAxes,
       ha="left" if left else "right", va="top", fontsize=6, color=INK2)
a.set_xlim(-0.6, 2.6)
panel_label(a, "a")
print(f"example: {site} {gene} {var} n={len(df)} q={row.qval:.2e}")

# ---------------- b: edSites vs sample size ----------------
b = ax["b"]
s = summary
b.scatter(s.samples, s.edsites, s=14, c=[col(v) for v in s.brain], edgecolors="white", linewidths=0.5, zorder=3)
fit = np.polyfit(s.samples, s.edsites, 1)
xx = np.array([0, s.samples.max() * 1.05])
b.plot(xx, np.polyval(fit, xx), color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2)
r = np.corrcoef(s.samples, s.edsites)[0, 1]
label = set(s.nlargest(3, "edsites").tissue) | set(s[s.brain].nlargest(1, "edsites").tissue) | set(s.nlargest(1, "samples").tissue)
smallest = s.nsmallest(1, "edsites").tissue.iloc[0]
SHORT = {"Skin_Sun_Exposed_Lower_leg": "Skin (sun-exposed)", "Muscle_Skeletal": "Muscle", "Nerve_Tibial": "Nerve",
         "Brain_Cerebellar_Hemisphere": "Cerebellar hemisphere", "Kidney_Cortex": "Kidney cortex"}
for _, rr in s[s.tissue.isin(label)].iterrows():
    right = rr.samples > 0.6 * s.samples.max()
    b.annotate(SHORT.get(rr.tissue, pretty(rr.tissue)), (rr.samples, rr.edsites),
               xytext=(-4, 3) if right else ((6, -1) if rr.tissue == smallest else (4, 3)), textcoords="offset points",
               ha="right" if right else "left", va="center" if rr.tissue == smallest else "baseline", fontsize=5.5, color=INK2)
b.set_xlabel("Sample size"); b.set_ylabel("Number of edSites")
b.set_xlim(0, None); b.set_ylim(0, None)
b.text(0.04, 0.97, f"{len(s)} tissues, r = {r:.2f}", transform=b.transAxes, va="top", fontsize=6, color=INK2)
b.legend(handles=group_legend, frameon=False, loc="lower right", fontsize=6, handletextpad=0.2)
b.set_title("edSites per tissue", color=INK)
panel_label(b, "b")

# ---------------- c: lead variant position relative to the editing site ----------------
cc = ax["c"]
bins = np.linspace(-100, 100, 81)
for br in [False, True]:
    d = lead.loc[lead.brain == br, "start_distance"] / 1000
    h, e = np.histogram(d[(d >= -100) & (d <= 100)], bins=bins)
    cc.plot((e[:-1] + e[1:]) / 2, h / h.sum(), color=col(br), lw=1.3, label="Brain" if br else "Other tissues")
cc.axvline(0, color=INK2, lw=0.6, ls=(0, (3, 2)))
cc.set_xlabel("Lead variant position relative to\nediting site (kb)"); cc.set_ylabel("Fraction of edSites")
cc.set_xlim(-100, 100)
within = (lead.start_distance.abs() <= 10000).mean()
cc.text(0.03, 0.97, f"{within:.0%} within\n±10 kb", transform=cc.transAxes, ha="left", va="top", fontsize=6, color=INK2)
cc.legend(frameon=False, loc="center right", fontsize=6, handlelength=1.2)
cc.set_title("Lead variants near edited sites", color=INK)
panel_label(cc, "c")

# ---------------- d: effect size vs minor allele frequency (all edSites) ----------------
dd = ax["d"]
maf = np.minimum(lead.af, 1 - lead.af)
eff = lead.slope.abs()
hb = dd.hexbin(maf, eff, gridsize=38, extent=(0.05, 0.5, 0, 2.5), bins="log", cmap=SEQ, mincnt=1, linewidths=0)
mb = np.linspace(0.05, 0.5, 10)
med = [eff[(maf >= lo) & (maf < hi)].median() for lo, hi in zip(mb[:-1], mb[1:])]
dd.plot((mb[:-1] + mb[1:]) / 2, med, color=INK, lw=1.4,
        path_effects=[matplotlib.patheffects.withStroke(linewidth=3, foreground="white")], label="Median")
dd.set_xlim(0.05, 0.5); dd.set_ylim(0, 2.5)
dd.set_xlabel("Minor allele frequency"); dd.set_ylabel("|Effect size| (normalized editing)")
cb = fig.colorbar(hb, ax=dd, fraction=0.05, pad=0.03)
cb.set_label("edSites", fontsize=6); cb.ax.tick_params(labelsize=5.5, width=0.5, length=2); cb.outline.set_linewidth(0.4)
dd.legend(frameon=False, loc="upper right", fontsize=6, handlelength=1.2)
dd.set_title("Effect size vs allele frequency", color=INK)
panel_label(dd, "d")

# ---------------- e: sharing of edQTL effects between tissues ----------------
ee = ax["e"]
tissues = list(summary.tissue)
R = pd.DataFrame(np.nan, index=tissues, columns=tissues)
for (dsc, tgt), grp in cross.groupby(["discovery", "target"]):
    if len(grp) >= 30:
        R.loc[dsc, tgt] = np.corrcoef(grp.slope_discovery, grp.slope_target)[0, 1]
S = (R + R.T) / 2                                                    # symmetric: average of both directions
S = S.fillna(R).fillna(R.T)
np.fill_diagonal(S.values, 1.0)
D = 1 - S.fillna(S.stack().median()).values
np.fill_diagonal(D, 0)
o = leaves_list(linkage(squareform(D, checks=False), "average"))
So = S.iloc[o, o]
vmin = np.floor(np.nanmin(So.values) * 10) / 10
im = ee.imshow(So.values, cmap=SEQ, vmin=vmin, vmax=1, interpolation="none", aspect="equal")
brain_o = np.array([tissues[i].startswith("Brain_") for i in o])
for k, br in enumerate(brain_o):                                     # group strip on top
    ee.add_patch(plt.Rectangle((k - 0.5, -2.6), 1, 1.6, color=col(br), clip_on=False, lw=0))
ee.set_xlim(-0.5, len(o) - 0.5); ee.set_ylim(len(o) - 0.5, -2.8)
ee.set_xticks([]); ee.set_yticks([])
for sp in ee.spines.values():
    sp.set_visible(False)
cb2 = fig.colorbar(im, ax=ee, fraction=0.05, pad=0.03)
cb2.set_label("Effect correlation (r)", fontsize=6); cb2.ax.tick_params(labelsize=5.5, width=0.5, length=2)
cb2.outline.set_linewidth(0.4)
conc = (np.sign(cross.slope_discovery) == np.sign(cross.slope_target)).mean()
ee.set_xlabel(f"{len(o)} tissues, clustered\n{conc:.0%} of effects in the same direction", labelpad=4)
ee.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)],
          labels=["Other tissues", "Brain"], frameon=False, fontsize=5.5, loc="upper center",
          bbox_to_anchor=(0.5, 1.13), ncol=2, handlelength=0.8, handletextpad=0.3, columnspacing=0.8)
ee.set_title("Sharing of edQTL effects", color=INK, pad=14)
panel_label(ee, "e", x=-0.1)
print(f"sharing: median pairwise r {np.nanmedian(S.values[np.triu_indices(len(S), 1)]):.2f}, same direction {conc:.3f}, "
      f"{len(cross):,} lookups")

# ---------------- f: allelic heterogeneity vs sample size ----------------
ff = ax["f"]
s["multi_frac"] = s.multi / s.edsites
ff.scatter(s.samples, s.multi_frac, s=14, c=[col(v) for v in s.brain], edgecolors="white", linewidths=0.5, zorder=3)
ff.set_xlabel("Sample size"); ff.set_ylabel("edSites with ≥2 edQTLs")
ff.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
ff.set_xlim(0, None); ff.set_ylim(0, None)
tot = s.multi.sum() / s.edsites.sum()
ff.text(0.04, 0.97, f"all tissues: {tot:.1%}", transform=ff.transAxes, va="top", fontsize=6, color=INK2)
ff.legend(handles=group_legend, frameon=False, loc="lower right", fontsize=6, handletextpad=0.2)
ff.set_title("Multiple independent edQTLs", color=INK)
panel_label(ff, "f")
print(f"multiple independent: {tot:.3f} overall")

for k in "abcf":
    ax[k].grid(axis="y", color=GRID, lw=0.4); ax[k].set_axisbelow(True)

fig.savefig(OUT, dpi=600)
print("saved", OUT)
