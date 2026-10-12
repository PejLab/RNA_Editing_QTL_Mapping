import os
import sys
import pandas as pd
from collections import defaultdict, Counter
configfile: 'config/config.yaml'


bam_extension = config["bam_extension"]
min_cov = config["min_coverage"]
min_sam = config["min_samples"]
maf = config.get('maf_threshold', 0.05)
geno_thresh = config.get('geno_threshold', 0.05)
hwe_thresh = config.get('hwe_threshold', 1e-6)
genotype_vcf = config.get('genotype_vcf', "data/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz")
prune_maf = config.get('prune_maf', 0.05)
prune_window = config.get('prune_window', 200)
prune_step = config.get('prune_step', 100)
prune_r2 = config.get('prune_r2', 0.1)
geno_pcs = config.get('geno_pcs', 3)
pheno_pcs = config.get('pheno_pcs', 5)
cis_window = config.get('cis_window', 100000)
fdr = config.get('fdr', 0.05)
plink_threads = config.get('plink_threads', 4)
normalization = config.get('normalization', 'quantile_int')
# Conda environments (used with --use-conda): a YAML file or the name of an existing environment
conda_env = config.get('conda_env', 'envs/environment.yml')
tensorqtl_env = config.get('tensorqtl_env', conda_env)

# Optional: regress host-gene expression out of editing levels (Li et al. 2022)
expression_covariate = config.get('expression_covariate', False)
gene_tpm = config.get('gene_tpm', '')
gene_annotation = config.get('gene_annotation', 'src/edit_site_gene_annotation.tsv.gz')
_bed_prefix = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt"
pheno_bed = _bed_prefix + (".qqnorm.expr_resid.bed.gz" if expression_covariate else ".qqnorm.bed.gz")   # used by covariates + tensorQTL


# Chrs
chrs = list(range(1, 23))  


# Reading annotation file
bam_anno = pd.read_csv(config["annotation_file"], sep='\t')


# Extracting tissue-bam pairs and only keep samples with BAMs present in bam dir
valid_pairs = [
    (tissue, sample)
    for tissue, sample in zip(bam_anno["tissue_cln"], bam_anno["sample_id"])
    if os.path.exists(f"data/bam/{sample}.{bam_extension}")
]


# Drop tissues with fewer samples than min_samples (no site can pass, so prep_phen would fail)
tissue_counts = Counter(t for t, s in valid_pairs)
small_tissues = sorted(t for t, n in tissue_counts.items() if n < min_sam)
if small_tissues:
    print(f"Skipping {len(small_tissues)} tissue(s) with < {min_sam} samples: {', '.join(small_tissues)}", file=sys.stderr)
valid_pairs = [(t, s) for t, s in valid_pairs if t not in small_tissues]


# Extract unique valid tissue and sort
tissues = sorted(set(t for t, s in valid_pairs))


# Making tissues > samples dic (for rule 2)
tissue_to_samples = defaultdict(list)
for tissue, sample in valid_pairs:
    tissue_to_samples[tissue].append(sample)


# Outputs - rule 1
query_editing_level_outputs = [
    f"intermediate/edit_level/{tissue}/{sample}.rnaeditlevel.gz"
    for tissue, sample in valid_pairs
]


# Outputs - rule 2
shared_sample_site_outputs = [
    f"intermediate/edit_matrix/{t}/{t}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt"
    for t in tissues
]


# Outputs - rule 3
prep_phen_outputs = [
    f"intermediate/edit_matrix/{t}/{t}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.phen_chr{n}"
    for t in tissues for n in range(1, 23)
]
prep_phen_qqnorm_outputs = [
    f"intermediate/edit_matrix/{t}/{t}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm_chr{n}"
    for t in tissues for n in range(1, 23)
]
prep_phen_misc_outputs = [
    f"intermediate/edit_matrix/{t}/{t}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt{suffix}"
    for t in tissues for suffix in [".ave", ".prepare.sh"]
]


# Outputs - rule 4
combine_to_bed_outputs = [
    f"intermediate/edit_matrix/{t}/{t}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.bed.gz"
    for t in tissues
]


# Outputs - rule 5
preparing_genotype_outputs = [
    f"data/genotype/{t}/{t}.geno.clean.maf{maf}.{ext}"
    for t in tissues
    for ext in ["bed", "bim", "fam"]
]
preparing_genotype_sample_list_outputs = [
    f"data/genotype/{t}/{t}.samples_to_keep.list"
    for t in tissues
]


# Outputs - rule 6
pruning_for_covar_outputs = [
    f"data/genotype/{t}/{t}.geno.pruned.vcf.gz" for t in tissues
]


# Outputs - rule 7
covariates_outputs = [
    f"data/genotype/{t}/{t}.covars.txt" for t in tissues
]
covariates_sex_age_outputs =[
    f"data/genotype/{t}/{t}.sex_age_info.txt" for t in tissues
]


# Outputs - rule 8
tensorqtl_cis_outputs = [
    f"results/{t}/{t}.cis_edQTLs.txt.gz" for t in tissues
]


# Outputs - rule 9
tensorqtl_cis_independent_outputs = [
    f"results/{t}/{t}.cis_independent_edQTLs.txt.gz" for t in tissues
]


# Outputs - rule 10
tensorqtl_cis_nominal_outputs = [
    f"results/{t}/nominal/{t}.cis_qtl_pairs.chr{n}.parquet"
    for t in tissues for n in range(1, 23)
]


# Outputs - rule 11
tensorqtl_all_signif_outputs = [
    f"results/{t}/{t}.cis_edqtl_signif.txt.gz" for t in tissues
]


# Outputs - rule 12
tensorqtl_all_cis_pvals_outputs = [
    f"results/{t}/{t}.cis_edqtl_all_pvals.txt.gz" for t in tissues
]


### END