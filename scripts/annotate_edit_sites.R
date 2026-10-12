#!/usr/bin/env Rscript
# Annotate editing sites with the genes they fall in, using an Ensembl GTF. As in Li et al. 2022, genes on both
# strands are kept (GTEx RNA-seq is not strand-specific; for cis-NAT sites both the sense and antisense gene count);
# column `sense` marks genes on the same strand as the edit.
# Site IDs match the pipeline: chrN_<1-based position> (BED end column).
#   Rscript scripts/annotate_edit_sites.R \
#       --edit_sites src/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed \
#       --gtf Homo_sapiens.GRCh38.110.chr.gtf \
#       --output src/edit_site_gene_annotation.tsv.gz

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(rtracklayer)
    library(GenomicRanges)
    library(argparse)
})

parser <- ArgumentParser(description = "Annotate editing sites with overlapping genes")
parser$add_argument("--edit_sites", required = TRUE, help = "BED file of editing sites (chrom, start, end, name, score, strand)")
parser$add_argument("--gtf", required = TRUE, help = "Ensembl GTF (chromosomes without 'chr' prefix)")
parser$add_argument("--output", required = TRUE, help = "Output TSV (.gz to compress)")
args <- parser$parse_args()

# Editing sites: BED is 0-based, so the 1-based position of the site is the end column
sites <- fread(args$edit_sites, header = FALSE) %>%
    select(chr_edsite = V1, pos_edsite = V3, strand_edsite = V6)

# Genes from the GTF
genes <- readGFF(args$gtf, filter = list(type = "gene"),
                 tags = c("gene_id", "gene_name", "gene_biotype")) %>%
    as.data.frame() %>%
    mutate(chr_gene = paste0("chr", seqid)) %>%
    select(chr_gene, start_gene = start, end_gene = end, strand_gene = strand, gene_id, gene_name, gene_biotype)

# Overlap each site (1 bp) with genes on either strand
hits <- findOverlaps(
    GRanges(sites$chr_edsite, IRanges(sites$pos_edsite, sites$pos_edsite), strand = sites$strand_edsite),
    GRanges(genes$chr_gene, IRanges(genes$start_gene, genes$end_gene), strand = genes$strand_gene),
    ignore.strand = TRUE)

annot <- cbind(sites[queryHits(hits), ], genes[subjectHits(hits), ]) %>%
    mutate(chr_pos_edsite = paste0(chr_edsite, "_", pos_edsite),
           sense = strand_edsite == strand_gene) %>%
    select(chr_pos_edsite, chr_edsite, pos_edsite, strand_edsite, gene_id, gene_name, gene_biotype, sense,
           chr_gene, start_gene, end_gene, strand_gene)

fwrite(annot, args$output, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("%d sites annotated with %d site-gene pairs (%d antisense)\n", uniqueN(annot$chr_pos_edsite), nrow(annot), sum(!annot$sense)))
