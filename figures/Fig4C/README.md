# Figure 4C. Bulk RNA-seq integrated with the Sun et al. single-cell MK atlas

Script: [`Fig4C_bulk_vs_sun_atlas.R`](Fig4C_bulk_vs_sun_atlas.R)

## What this analysis does

Sorted MHC II+ and MHC II- bone marrow MKs were profiled by bulk RNA-seq (n = 4 per group). This analysis relates those profiles to the four MK subsets defined by Sun et al. (2021): platelet-producing, HSC-niche, inflammatory/immune, and actively cycling.

The comparison is quantitative. Each reference subset is collapsed to pseudobulk, each bulk replicate is assigned to the subset it correlates with most strongly, and the bulk MHC II signatures are scored on every reference cell. A positive control checks that the mapping is specific.

## Steps

1. **Annotate the reference.** The four clusters at `RNA_snn_res.0.1` are mapped to the four subsets in `SUBSET_MAP`. The marker evidence for the map is written to `qc_cluster_marker_means.csv` on every run.
2. **Pseudobulk.** Raw counts are summed within each subset. Cells in each subset are split at random into three pseudo-replicates to capture within-subset variance.
3. **Unified gene panel.** The top 12 markers per subset (`FindAllMarkers`, Wilcoxon, adjusted P < 0.05) are combined with 15 MHC class II and co-stimulatory genes.
4. **Panel heatmaps.** The panel is shown as row z-scored log2 CPM (edgeR) for the reference pseudobulks and for the bulk samples.
5. **Mapping.** Each bulk replicate is correlated with each subset pseudobulk (Spearman, 2,000 shared highly variable genes) and assigned to its best match.
6. **Signature scoring.** The bulk MHC II+ and MHC II- DE signatures are scored on every reference cell with `AddModuleScore`.
7. **Validation control.** The subset pseudo-replicates are passed back through step 5. A specific mapping returns each one to its own subset.

## Inputs

Three files in `data/`: `mouse_mk.rds`, `counts_matrix.tsv` and `DE_genes_sig.tsv`. See [`data/README.md`](../../data/README.md).

## Outputs

Written to `results/Fig4C/`.

| File | Content |
|---|---|
| `pseudobulk_panel_heatmap.pdf` | Unified panel across the reference subset pseudobulks |
| `bulk_panel_heatmap.pdf` | Unified panel across the bulk MHC II+ and MHC II- samples |
| `bulk_vs_subset_correlation.pdf` | Spearman correlation of each bulk replicate to each subset |
| `signature_scores_violin.pdf` | MHC II+ and MHC II- signature scores by subset |
| `signature_scores_umap.pdf` | The same scores on the reference UMAP |
| `mapping_validation_control.pdf` | Validation control |
| `qc_subset_marker_dotplot.pdf`, `qc_cluster_marker_means.csv` | Evidence for the subset annotation |
| `bulk_vs_subset_spearman.csv`, `bulk_nearest_subset.csv` | Mapping results |
| `mapping_validation_control.csv` | Validation control results |
| `sun_pseudobulk_counts.csv`, `sun_subset_markers.csv`, `unified_panel_genes.csv` | Intermediate tables |
| `sessionInfo.txt` | R environment for the run |

Set `EXPORT_TIFF <- TRUE` in the config block to also write 300 dpi TIFFs.

## Running

From the repository root:

```bash
Rscript figures/Fig4C/Fig4C_bulk_vs_sun_atlas.R
```

Input and output locations can be set without editing the script, through the environment variables `MK_SC_RDS`, `MK_BULK_COUNTS`, `MK_DE_TABLE` and `MK_OUTDIR`.

## Requirements

| Package | Version tested |
|---|---|
| R | 4.6.0 |
| Seurat | 5.5.1 |
| SeuratObject | 5.4.0 |
| edgeR | 4.10.1 |
| pheatmap | 1.0.13 |
| ggplot2 | 4.0.3 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| RColorBrewer | 1.1 |

```r
install.packages(c("Seurat", "pheatmap", "RColorBrewer", "ggplot2", "dplyr", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("edgeR")
```

The full environment is recorded in [`environment/sessionInfo_Fig4C.txt`](../../environment/sessionInfo_Fig4C.txt).

## Reproducibility

The seed is fixed (`SEED <- 13579`), so the pseudo-replicate split and the module scores are identical between runs. The subset annotation is fixed by cluster ID. If the reference object is rebuilt and the cluster IDs change, check `qc_cluster_marker_means.csv` and update `SUBSET_MAP` before using any downstream output.
