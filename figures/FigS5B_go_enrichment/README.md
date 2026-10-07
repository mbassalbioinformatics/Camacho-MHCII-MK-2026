# Fig S5B. GO biological process enrichment

Script: [`FigS5B_go_enrichment.R`](FigS5B_go_enrichment.R)

## What this analysis does

The genes higher in MHC II+ MKs and the genes higher in MHC II- MKs (DESeq2, absolute log2 fold change > 1, adjusted P < 0.05) are tested separately for enrichment of Gene Ontology biological processes.

## Steps

1. The DESeq2 model is rebuilt with [`R/bulk_rnaseq_helpers.R`](../../R/bulk_rnaseq_helpers.R), identical to Fig 4A.
2. Gene symbols are mapped to Entrez IDs with `bitr` (org.Mm.eg.db).
3. Each gene set is tested with `clusterProfiler::enrichGO`: ontology BP, Benjamini-Hochberg correction, P and q < 0.05. The background is the clusterProfiler default (all annotated mouse genes).
4. **Fig S5B.** The top 20 terms per group are shown as dot plots side by side.

## Outputs

Written to `results/FigS5B_go_enrichment/`.

| File | Content |
|---|---|
| `FigS5B_GO_BP_dotplot.pdf` | Fig S5B |
| `GO_BP_MHCpos.tsv`, `GO_BP_MHCneg.tsv` | Full enrichment results for each group |
| `sessionInfo.txt` | R environment for the run |

## Reproducibility

The input gene lists are identical to the published analysis. GO annotations, however, are updated with each Bioconductor release. With the versions tested here, the top terms are similar but not identical to the published panel, which was generated with an earlier annotation release. The themes reported in the paper are all recovered: phagocytosis, interferon responses and cell adhesion for MHC II+ genes, and leukocyte migration, chemotaxis and coagulation for MHC II- genes.

## Running

```bash
Rscript figures/FigS5B_go_enrichment/FigS5B_go_enrichment.R
```

Input: `data/counts_matrix.tsv` (see [`data/README.md`](../../data/README.md)).

## Requirements

Tested with R 4.6.0, DESeq2 1.52.0, clusterProfiler 4.20.0, org.Mm.eg.db 3.23.0, GO.db 3.23.1 and enrichplot 1.32.1 (Bioconductor 3.23). The full environment is in [`environment/sessionInfo_FigS5B.txt`](../../environment/sessionInfo_FigS5B.txt).
