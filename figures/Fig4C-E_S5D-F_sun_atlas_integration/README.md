# Fig 4C-E and Fig S5D-F. Bulk RNA-seq integrated with the Sun et al. single-cell MK atlas

Script: [`Fig4C-E_S5D-F_sun_atlas_integration.R`](Fig4C-E_S5D-F_sun_atlas_integration.R)

## What this analysis does

This analysis relates the MHC II+ and MHC II- bulk profiles to the four MK subsets defined by Sun et al. (2021): platelet-producing, HSC-niche, inflammatory/immune, and actively cycling.

The comparison is quantitative. Each reference subset is collapsed to pseudobulk, each bulk replicate is assigned to the subset it correlates with most strongly, and the bulk MHC II signatures are scored on every reference cell. A positive control checks that the mapping is specific.

## Steps

1. **Annotate the reference.** The four clusters at `RNA_snn_res.0.1` are mapped to the four subsets in `SUBSET_MAP`. The marker evidence for the map is written to `qc_cluster_marker_means.csv` on every run.
2. **Pseudobulk.** Raw counts are summed within each subset. Cells in each subset are split at random into three pseudo-replicates to capture within-subset variance.
3. **Unified gene panel.** The top 12 markers per subset (`FindAllMarkers`, Wilcoxon, adjusted P < 0.05) are combined with 15 MHC class II and co-stimulatory genes, giving 56 genes.
4. **Panel heatmaps (Fig 4C, 4D).** The panel is shown as row z-scored log2 CPM (edgeR) for the reference pseudobulks and for the bulk samples.
5. **Mapping (Fig 4E).** Each bulk replicate is correlated with each subset pseudobulk (Spearman, 2,000 shared highly variable genes) and assigned to its best match.
6. **Validation control (Fig S5D).** The subset pseudo-replicates are passed back through step 5. A specific mapping returns each one to its own subset (12 of 12 in the published analysis).
7. **Signature scoring (Fig S5E, S5F).** The bulk MHC II+ and MHC II- DE signatures are scored on every reference cell with `AddModuleScore`.

## Inputs

- `data/mouse_mk.rds` and `data/counts_matrix.tsv` (see [`data/README.md`](../../data/README.md)).
- `results/Fig4A_S5A_differential_expression/DE_genes_sig.tsv`, written by the [Fig 4A/S5A script](../Fig4A_S5A_differential_expression). Run that script first, or set `MK_DE_TABLE` to an existing copy.

## Outputs

Written to `results/Fig4C-E_S5D-F_sun_atlas_integration/`.

| File | Content |
|---|---|
| `Fig4C_pseudobulk_panel_heatmap.pdf` | Fig 4C |
| `Fig4D_bulk_panel_heatmap.pdf` | Fig 4D |
| `Fig4E_bulk_vs_subset_correlation.pdf` | Fig 4E |
| `FigS5D_mapping_validation_control.pdf` | Fig S5D |
| `FigS5E_signature_scores_violin.pdf` | Fig S5E |
| `FigS5F_signature_scores_umap.pdf` | Fig S5F |
| `qc_subset_marker_dotplot.pdf`, `qc_cluster_marker_means.csv` | Evidence for the subset annotation |
| `bulk_vs_subset_spearman.csv`, `bulk_nearest_subset.csv` | Mapping results |
| `mapping_validation_control.csv` | Validation control results |
| `sun_pseudobulk_counts.csv`, `sun_subset_markers.csv`, `unified_panel_genes.csv` | Intermediate tables |
| `sessionInfo.txt` | R environment for the run |

Set `EXPORT_TIFF <- TRUE` in the config block to also write 300 dpi TIFFs.

## Running

From the repository root, after the Fig 4A/S5A script:

```bash
Rscript figures/Fig4C-E_S5D-F_sun_atlas_integration/Fig4C-E_S5D-F_sun_atlas_integration.R
```

Input and output locations can be set without editing the script, through the environment variables `MK_SC_RDS`, `MK_BULK_COUNTS`, `MK_DE_TABLE` and `MK_OUTDIR`.

## Requirements

Tested with the following versions. The full environment is in [`environment/sessionInfo_Fig4C-E_S5D-F.txt`](../../environment/sessionInfo_Fig4C-E_S5D-F.txt).

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

## Reproducibility

This script reproduces the published correlation, nearest-subset and validation-control tables exactly. The seed is fixed (`SEED <- 13579`), so the pseudo-replicate split and the module scores are identical between runs. The subset annotation is fixed by cluster ID. If the reference object is rebuilt and the cluster IDs change, check `qc_cluster_marker_means.csv` and update `SUBSET_MAP` before using any downstream output.
