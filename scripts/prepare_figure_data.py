#!/usr/bin/env python3
# Collect the summary tables used by make_summary_figure.py from the pipeline results (all tissues).
# Run once after the pipeline (heavy step: reads every tissue's nominal statistics); writes results/figures/data/.
#   python3 scripts/prepare_figure_data.py [n_processes]

import glob
import os
import sys
from multiprocessing import Pool
import pandas as pd
import pyarrow.parquet as pq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = f"{ROOT}/results"
OUT = f"{ROOT}/results/figures/data"
FDR = 0.05
os.makedirs(OUT, exist_ok=True)
TISSUES = sorted(os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{RES}/*/*.cis_edQTLs.txt.gz"))   # tissue result folders only


def tissue_tables(t):
    cis = pd.read_csv(f"{RES}/{t}/{t}.cis_edQTLs.txt.gz", sep="\t")
    ind = pd.read_csv(f"{RES}/{t}/{t}.cis_independent_edQTLs.txt.gz", sep="\t", usecols=["phenotype_id"])
    sig = cis[cis.qval <= FDR]
    n_ind = ind.groupby("phenotype_id").size().reindex(sig.phenotype_id).fillna(1)
    n_samples = len(glob.glob(f"{ROOT}/intermediate/edit_level/{t}/*.rnaeditlevel.gz"))
    summary = dict(tissue=t, samples=n_samples, sites_tested=len(cis), edsites=len(sig),
                   multi=int((n_ind >= 2).sum()))
    lead = sig[["phenotype_id", "variant_id", "start_distance", "af", "slope", "qval", "pval_beta"]].assign(tissue=t)
    return summary, lead


def lookup(args):
    """Effects in tissue `t` of every tissue's lead (site, variant) pairs, from t's nominal statistics."""
    t, pairs = args
    hits = []
    for f in sorted(glob.glob(f"{RES}/{t}/nominal/*.parquet")):
        d = pq.read_table(f, columns=["phenotype_id", "variant_id", "slope"]).to_pandas()
        hits.append(d.merge(pairs, on=["phenotype_id", "variant_id"], how="inner"))
    h = pd.concat(hits).rename(columns={"slope": "slope_target"}).assign(target=t)
    print(f"  {t}: {len(h):,} lead pairs found", flush=True)
    return h


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    summaries, leads = zip(*[tissue_tables(t) for t in TISSUES])
    summary = pd.DataFrame(summaries)
    lead = pd.concat(leads)
    summary.to_csv(f"{OUT}/tissue_summary.tsv", sep="\t", index=False)
    lead.to_csv(f"{OUT}/lead_edqtls.tsv.gz", sep="\t", index=False)
    print(f"{len(TISSUES)} tissues, {len(lead):,} edSites", flush=True)

    pairs = lead[["phenotype_id", "variant_id", "tissue", "slope"]].rename(columns={"tissue": "discovery", "slope": "slope_discovery"})
    with Pool(nproc) as pool:
        cross = pd.concat(pool.map(lookup, [(t, pairs) for t in TISSUES]))
    cross = cross[cross.discovery != cross.target]
    cross.to_csv(f"{OUT}/cross_tissue_effects.tsv.gz", sep="\t", index=False)
    print(f"cross-tissue lookups: {len(cross):,}", flush=True)
