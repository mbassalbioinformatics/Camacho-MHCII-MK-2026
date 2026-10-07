# Fig 4A and Fig S5A. Differential expression, MHC II+ vs MHC II- MKs

Script: [`Fig4A_S5A_differential_expression.R`](Fig4A_S5A_differential_expression.R)

## What this analysis does

Bulk RNA-seq of sorted MHC II+ and MHC II- bone marrow MKs (four pooled libraries per group) is compared with DESeq2. Genes with an absolute log2 fold change > 1 and an adjusted P < 0.05 are called differentially expressed.

## Steps

1. Gene symbols are taken from the `ENSMUSG..._Symbol` row names, immunoglobulin variable and alpha genes are removed (`^Ighv`, `^Igha`, `^Igkv`), duplicate symbols are summed, and counts are converted to integers. This is shared with the other bulk analyses through [`R/bulk_rnaseq_helpers.R`](../../R/bulk_rnaseq_helpers.R).
2. DESeq2 is fit with `~ Group`, MHC II- as the reference level, using the Wald test. Log2 fold changes are not shrunk.
3. **Fig 4A.** The volcano plot labels the ten most significant genes in each direction.
4. **Fig S5A.** The heatmap shows the 50 genes with the lowest adjusted P, as rlog values scaled by row.

## Outputs

Written to `results/Fig4A_S5A_differential_expression/`.

| File | Content |
|---|---|
| `Fig4A_volcano.pdf` | Fig 4A |
| `FigS5A_top50_DE_heatmap.pdf` | Fig S5A |
| `DE_all_genes.tsv` | DESeq2 results for every gene |
| `DE_genes_sig.tsv` | Significant genes with their direction (`Up_MHCpos` or `Up_MHCneg`); input to Fig 4C-E and S5D-F |
| `sessionInfo.txt` | R environment for the run |

Published result: 518 genes higher in MHC II+ and 235 higher in MHC II-.

## Running

```bash
Rscript figures/Fig4A_S5A_differential_expression/Fig4A_S5A_differential_expression.R
```

Input: `data/counts_matrix.tsv` (see [`data/README.md`](../../data/README.md)).

## Requirements

Tested with R 4.6.0, DESeq2 1.52.0, ggplot2 4.0.3, ggrepel 0.9.8 and pheatmap 1.0.13 (Bioconductor 3.23). The full environment is in [`environment/sessionInfo_Fig4A_S5A.txt`](../../environment/sessionInfo_Fig4A_S5A.txt).

## Reproducibility

This script reproduces the published DE table exactly: the same 753 genes, the same direction calls, and identical statistics to floating-point precision.
