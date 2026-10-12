#!/usr/bin/env python3
# Remove host-gene expression effects from editing levels, as in Li et al. 2022: for each editing site, regress
# its normalized editing level on the expression (log2 TPM + 1) of the gene(s) it lies in, and keep the residuals.
# Sites in overlapping genes (e.g. cis-NATs) use all of them; sites outside expressed genes are left unchanged.
# Residuals are re-normalized with the pipeline's normalization method.

import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm

parser = argparse.ArgumentParser(description="Regress host-gene expression out of editing levels")
parser.add_argument("--bed", required=True, help="Normalized phenotype BED (#Chr, start, end, ID, subjects...)")
parser.add_argument("--tpm", required=True, help="GTEx gene TPM table: .ft (feather) or .gct/.gct.gz")
parser.add_argument("--annotation", required=True, help="Sample annotation (sample_id, subject_id, tissue_cln)")
parser.add_argument("--tissue", required=True, help="Tissue name as in tissue_cln")
parser.add_argument("--gene_annotation", required=True, help="Site-to-gene table (chr_pos_edsite, gene_id)")
parser.add_argument("--norm", default="quantile_int", choices=["quantile_int", "zscore_qqnorm"],
                    help="Re-normalization of residuals: rank-based inverse normal (quantile_int) or z-score")
parser.add_argument("--output", required=True, help="Output BED (uncompressed)")
args = parser.parse_args()

bed = pd.read_csv(args.bed, sep="\t")
subjects = list(bed.columns[4:])

# Samples of this tissue, mapped to the subject IDs used in the BED
anno = pd.read_csv(args.annotation, sep="\t", dtype=str)
for col in ["sample_id", "subject_id", "tissue_cln"]:
    if col not in anno.columns:
        sys.exit(f"Error: {args.annotation} is missing column {col}")
anno = anno[(anno.tissue_cln == args.tissue) & anno.subject_id.isin(subjects)].drop_duplicates("subject_id")
sample_to_subject = dict(zip(anno.sample_id, anno.subject_id))

# Expression of those samples
if args.tpm.endswith(".ft"):
    import pyarrow.feather as feather
    table = feather.read_table(args.tpm, memory_map=True)          # memory-mapped: only used columns are loaded
    use = [s for s in sample_to_subject if s in set(table.column_names)]
    tpm = table.select(["Name"] + use).to_pandas()
else:
    header = pd.read_csv(args.tpm, sep="\t", skiprows=2, nrows=0).columns
    use = [s for s in sample_to_subject if s in header]
    tpm = pd.read_csv(args.tpm, sep="\t", skiprows=2, usecols=["Name"] + use)
tpm = tpm[~tpm.Name.str.endswith("_PAR_Y")]                       # duplicate chrY copies of PAR genes
tpm["gene_id"] = tpm.Name.str.split(".").str[0]                    # drop GENCODE version suffix
tpm = tpm.drop(columns="Name").set_index("gene_id")
tpm.columns = [sample_to_subject[s] for s in tpm.columns]
missing = [s for s in subjects if s not in tpm.columns]
if missing:
    sys.exit(f"Error: no expression for {len(missing)} subject(s) in {args.tissue}: {', '.join(missing[:10])}")
expr = np.log2(tpm[subjects] + 1)
expr = expr[expr.std(axis=1) > 0]                                  # expressed, variable genes only
print(f"{args.tissue}: {len(subjects)} subjects, expression for {len(expr):,} genes", flush=True)

# Host genes of each site
genes = pd.read_csv(args.gene_annotation, sep="\t", usecols=["chr_pos_edsite", "gene_id"])
genes = genes[genes.gene_id.isin(expr.index)]
site_genes = genes.groupby("chr_pos_edsite").gene_id.apply(lambda g: sorted(set(g)))


def renormalize(v):
    if args.norm == "quantile_int":
        return norm.ppf(rankdata(v) / (len(v) + 1))
    return (v - v.mean()) / v.std()


values = bed[subjects].to_numpy(float)
E = expr.to_numpy()
gene_row = {g: i for i, g in enumerate(expr.index)}
n_genes, r2 = [], []
for i, site in enumerate(bed.iloc[:, 3]):
    g = site_genes.get(site, [])
    n_genes.append(len(g))
    if not g:
        continue
    y = values[i]
    X = np.column_stack([np.ones(len(y))] + [E[gene_row[x]] for x in g])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    r2.append(1 - resid.var() / y.var())
    values[i] = renormalize(resid)

bed[subjects] = values
bed.to_csv(args.output, sep="\t", index=False, float_format="%g")
n_genes = np.array(n_genes)
print(f"sites: {len(n_genes):,} | with 1 host gene: {(n_genes == 1).sum():,} | with 2+: {(n_genes >= 2).sum():,} | "
      f"no expressed host gene (unchanged): {(n_genes == 0).sum():,}", flush=True)
if r2:
    print(f"variance explained by host-gene expression: median {np.median(r2):.3f}, >10% for {np.mean(np.array(r2) > 0.1):.1%} of sites",
          flush=True)
