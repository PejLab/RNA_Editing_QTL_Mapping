configfile: 'config/config.yaml'
include:"rules/globals.smk"


rule all:
    input:
        query_editing_level_outputs,
        shared_sample_site_outputs,
        prep_phen_outputs + prep_phen_qqnorm_outputs + prep_phen_misc_outputs,
        combine_to_bed_outputs,
        preparing_genotype_outputs + preparing_genotype_sample_list_outputs,
        pruning_for_covar_outputs,
        covariates_outputs + covariates_sex_age_outputs,
        tensorqtl_cis_outputs,
        tensorqtl_cis_independent_outputs,
        tensorqtl_cis_nominal_outputs,
        tensorqtl_all_signif_outputs,
        tensorqtl_all_cis_pvals_outputs


rule index_ref_genome:
    """Index the reference once, so parallel mpileup jobs don't all try to create the .fai at the same time."""
    input:
        config["ref_genome_fasta"]
    output:
        config["ref_genome_fasta"] + ".fai"
    conda:
        conda_env
    shell:
        "samtools faidx {input}"


rule query_editing_level:
    input:
        edit_ref = config["edit_sites_bed"],
        ref_genome = config["ref_genome_fasta"],
        ref_genome_fai = config["ref_genome_fasta"] + ".fai",
        bam = "data/bam/{sample}." + bam_extension
    output:
        out = "intermediate/edit_level/{tissue}/{sample}.rnaeditlevel.gz"
    log:
        "logs/query_editing_level/{tissue}_{sample}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p intermediate/edit_level/{wildcards.tissue} logs/query_editing_level/
        perl scripts/query_editing_level.pl \
        --edit_sites {input.edit_ref} \
        --ref {input.ref_genome} \
        --bam {input.bam} \
        --output {output.out} &> {log}
        """


rule shared_sample_site:
    input:
        lambda wildcards: [
            f"intermediate/edit_level/{wildcards.tissue}/{sample}.rnaeditlevel.gz"
            for sample in tissue_to_samples[wildcards.tissue]
        ]
    output:
        out = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt"
    params:
        min_cov = min_cov,
        min_sam = min_sam
    log:
        "logs/shared_sample_site/{tissue}.shared_sample_site.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p intermediate/edit_matrix/{wildcards.tissue} logs/shared_sample_site/
        perl scripts/shared_samples_sites_matrix.pl \
            --path_to_edit_files intermediate/edit_level/{wildcards.tissue}/ \
            --output_dir intermediate/edit_matrix/{wildcards.tissue}/ \
            --tissue {wildcards.tissue} \
            --min_coverage {params.min_cov} \
            --min_samples {params.min_sam} \
            &> {log}
        """


rule prep_phen:
    input:
        matrix = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt"
    output:
        phen_chr = [
            f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.phen_chr{chrn}"
            for chrn in chrs
        ],
        qqnorm_chr = [
            f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm_chr{chrn}"
            for chrn in chrs
        ],
        ave = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.ave",
        script = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.prepare.sh"
    params:
        norm = normalization
    log:
        f"logs/prep_phen/{{tissue}}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p logs/prep_phen/
        python3 scripts/prepare_phenotype_table_py3.py -p 0 -n {params.norm} {input.matrix} &> {log} # -p 0: skip the PCA (.PCs file), covariate PCs are made in the covariates rule
        """


rule combine_to_bed:
    input:
        qqnorm_chr = [
            f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm_chr{chrn}"
            for chrn in chrs
        ],
    output:
        bed = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.bed.gz",
        index = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.bed.gz.tbi"
    params:
        input_dir = lambda wildcards: f"intermediate/edit_matrix/{wildcards.tissue}",
        tissue = f"{{tissue}}",
        min_cov = min_cov,
        min_sam = min_sam
    log:
        "logs/combine_to_bed/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p logs/combine_to_bed/
        python3 scripts/combine_to_bed.py \
                --input_dir {params.input_dir} \
                --tissue {params.tissue} \
                --min_coverage {params.min_cov} \
                --min_samples {params.min_sam} &> {log}
        """


rule regress_expression:
    """Regress host-gene expression out of each editing site's levels (optional: expression_covariate)."""
    input:
        bed = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.bed.gz",
        tpm = gene_tpm,
        annotation = config["annotation_file"],
        genes = gene_annotation
    output:
        bed = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.expr_resid.bed.gz",
        index = f"intermediate/edit_matrix/{{tissue}}/{{tissue}}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt.qqnorm.expr_resid.bed.gz.tbi"
    params:
        bed = lambda wildcards, output: output.bed[:-3],
        norm = normalization
    log:
        "logs/regress_expression/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        python3 scripts/regress_expression.py \
            --bed {input.bed} \
            --tpm {input.tpm} \
            --annotation {input.annotation} \
            --tissue {wildcards.tissue} \
            --gene_annotation {input.genes} \
            --norm {params.norm} \
            --output {params.bed} &> {log}
        bgzip -f {params.bed}
        tabix -f -p bed {output.bed}
        """


rule preparing_genotype:
    input:
        vcf = genotype_vcf,
        ref_genome = config["ref_genome_fasta"],
        annotation_file = config["annotation_file"]
    output:
        multiext("data/genotype/{tissue}/{tissue}.geno.clean.maf" + str(maf), ".bed", ".bim", ".fam"),
        sample_list = "data/genotype/{tissue}/{tissue}.samples_to_keep.list"
    params:
        maf = maf,
        geno_thresh = geno_thresh,
        hwe_thresh = hwe_thresh,
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    threads: plink_threads
    log:
        "logs/preparing_genotype/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p $(dirname {output[0]}) logs/preparing_genotype/

        python3 scripts/subset_annotation.py \
            --annotation {input.annotation_file} \
            --tissue {wildcards.tissue} \
            --output {output.sample_list}

        plink2 \
            --vcf {input.vcf} \
            --fa {input.ref_genome} \
            --ref-from-fa force \
            --keep {output.sample_list} \
            --nonfounders \
            --chr 1-22 \
            --snps-only 'just-acgt' \
            --maf {params.maf} \
            --geno {params.geno_thresh} \
            --hwe {params.hwe_thresh} \
            --make-bed \
            --output-chr chrM \
            --threads {threads} \
            --out {params.prefix} &> {log}
        """


rule pruning_for_covar:
    input:
        multiext("data/genotype/{tissue}/{tissue}.geno.clean.maf" + str(maf), ".bed", ".bim", ".fam")
    output:
        prunned_vcf = "data/genotype/{tissue}/{tissue}.geno.pruned.vcf.gz"
    params:
        geno_prefix = lambda wildcards: f"data/genotype/{wildcards.tissue}/{wildcards.tissue}.geno.clean.maf{maf}",
        pruned_prefix = lambda wildcards: f"data/genotype/{wildcards.tissue}/{wildcards.tissue}.geno.pruning_list",
        prune_maf = prune_maf,
        prune_window = prune_window,
        prune_step = prune_step,
        prune_r2 = prune_r2
    threads: plink_threads
    log:
        "logs/pruning_for_covar/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p logs/pruning_for_covar/

        plink2 \
            --bfile {params.geno_prefix} \
            --maf {params.prune_maf} \
            --indep-pairwise {params.prune_window} {params.prune_step} {params.prune_r2} \
            --bad-ld \
            --threads {threads} \
            --out {params.pruned_prefix} &> {log}

        plink2 \
            --bfile {params.geno_prefix} \
            --extract {params.pruned_prefix}.prune.in \
            --export vcf bgz id-paste=iid \
            --threads {threads} \
            --out {params.pruned_prefix} &>> {log}

        mv {params.pruned_prefix}.vcf.gz {output.prunned_vcf}
        """


rule covariates:
    input:
        vcf_pruned = f"data/genotype/{{tissue}}/{{tissue}}.geno.pruned.vcf.gz",
        bed = pheno_bed,
        annotation_file = config["annotation_file"],
        subject_covariates = config["subject_covariates"]
    output:
        out = f"data/genotype/{{tissue}}/{{tissue}}.covars.txt",
        sex_age = f"data/genotype/{{tissue}}/{{tissue}}.sex_age_info.txt"
    params:
        geno_pcs = geno_pcs,
        pheno_pcs = pheno_pcs
    log:
        "logs/covariates/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        python3 scripts/subset_annotation.py \
            --annotation {input.annotation_file} \
            --covariates {input.subject_covariates} \
            --tissue {wildcards.tissue} \
            --output {output.sex_age}
        Rscript scripts/covariates.R \
                --vcf {input.vcf_pruned} \
                --bed {input.bed} \
                --sex_age_info {output.sex_age} \
                --geno_pcs {params.geno_pcs} \
                --pheno_pcs {params.pheno_pcs} \
                --out {output.out} &> {log}
        """


rule tensorqtl_cis:
    input:
        geno = multiext("data/genotype/{tissue}/{tissue}.geno.clean.maf" + str(maf), ".bed", ".bim", ".fam"),
        bed = pheno_bed,
        bedi = pheno_bed + ".tbi",
        covar = f"data/genotype/{{tissue}}/{{tissue}}.covars.txt",
    output:
        out = f"results/{{tissue}}/{{tissue}}.cis_edQTLs.txt.gz"
    params:
        geno_prefix = lambda wildcards: f"data/genotype/{wildcards.tissue}/{wildcards.tissue}.geno.clean.maf{maf}",
        window = cis_window
    resources:
        tensorqtl_slots = 1
    log:
        "logs/tensorqtl_cis/{tissue}.log"
    conda:
        tensorqtl_env
    shell:
        """
        mkdir -p results/{wildcards.tissue} logs/tensorqtl_cis/
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output.out} \
            --covariates {input.covar} \
            --mode cis \
            --window {params.window} &> {log}
        """


rule tensorqtl_cis_independent:
    input:
        geno = multiext("data/genotype/{tissue}/{tissue}.geno.clean.maf" + str(maf), ".bed", ".bim", ".fam"),
        bed = pheno_bed,
        bedi = pheno_bed + ".tbi",
        covar = f"data/genotype/{{tissue}}/{{tissue}}.covars.txt",
        cis_out = f"results/{{tissue}}/{{tissue}}.cis_edQTLs.txt.gz"
    output:
        out = f"results/{{tissue}}/{{tissue}}.cis_independent_edQTLs.txt.gz"
    params:
        geno_prefix = lambda wildcards: f"data/genotype/{wildcards.tissue}/{wildcards.tissue}.geno.clean.maf{maf}",
        window = cis_window
    resources:
        tensorqtl_slots = 1
    log:
        "logs/tensorqtl_cis_independent/{tissue}.log"
    conda:
        tensorqtl_env
    shell:
        """
        mkdir -p logs/tensorqtl_cis_independent/
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output.out} \
            --covariates {input.covar} \
            --cis_output {input.cis_out} \
            --mode cis_independent \
            --window {params.window} &> {log}
        """


rule tensorqtl_cis_nominal:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("data/genotype/{tissue}/{tissue}.geno.clean.maf" + str(maf), ".bed", ".bim", ".fam"),
        bed = pheno_bed,
        bedi = pheno_bed + ".tbi",
        covar = f"data/genotype/{{tissue}}/{{tissue}}.covars.txt"
    output:
        expand("results/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=range(1, 23))
    params:
        geno_prefix = lambda wildcards: f"data/genotype/{wildcards.tissue}/{wildcards.tissue}.geno.clean.maf{maf}",
        outdir = "results/{tissue}/nominal",
        out_prefix = "{tissue}",
        window = cis_window
    resources:
        tensorqtl_slots = 1
    log:
        "logs/tensorqtl_cis_nominal/{tissue}.log"
    conda:
        tensorqtl_env
    shell:
        """
        mkdir -p {params.outdir}
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {params.out_prefix} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --mode cis_nominal \
            --window {params.window} &> {log}
        """


rule tensorqtl_all_signif:
    """Extract all significant cis SNP-gene pairs."""
    input:
        perm = "results/{tissue}/{tissue}.cis_edQTLs.txt.gz",
        nom = expand("results/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=chrs)
        # nom = [
        #     f"results/{tissue}/nominal/{tissue}.cis_qtl_pairs.chr{chrn}.parquet"
        #     for chrn in chrs
        #     ]
    output:
        out = "results/{tissue}/{tissue}.cis_edqtl_signif.txt.gz"
    params:
        nom_prefix = "results/{tissue}/nominal/{tissue}",
        fdr = fdr
    log:
        "logs/tensorqtl_all_signif/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p logs/tensorqtl_all_signif
        python3 scripts/tensorqtl_all_signif.py \
            {input.perm} \
            {params.nom_prefix} \
            {output.out} \
            --fdr {params.fdr} &> {log}
        """


rule tensorqtl_all_cis_pvals:
    """Extract p-values for all tested cis-window SNPs per gene."""
    input:
        nom = expand("results/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=chrs)
        # nom = [
        #     f"results/{tissue}/nominal/{tissue}.cis_qtl_pairs.chr{chrn}.parquet"
        #     for chrn in chrs
        #     ]
    output:
        out = "results/{tissue}/{tissue}.cis_edqtl_all_pvals.txt.gz"
    params:
        nom_dir = "results/{tissue}/nominal"
    log:
        "logs/tensorqtl_all_cis_pvals/{tissue}.log"
    conda:
        conda_env
    shell:
        """
        mkdir -p logs/tensorqtl_all_cis_pvals
        python3 scripts/tensorqtl_all_cis_pvals.py \
            {params.nom_dir} \
            {output.out} &> {log}
        """


### END