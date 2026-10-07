# Input data

Input files go in this folder. None of them are tracked in git.

| File | Used by | Content | Source |
|---|---|---|---|
| `counts_matrix.tsv` | All bulk RNA-seq analyses | RSEM expected counts, genes by 8 samples, tab-separated. Row names have the form `ENSMUSG..._Symbol`. Columns 1 to 4 are MHC II+ and columns 5 to 8 are MHC II-. | Generated from the bulk RNA-seq deposited at GEO (GSE328415) |
| `mouse_mk.rds` | Fig 4C-E, S5D-F | Seurat object of the Sun et al. (2021) bone marrow MK scRNA-seq. It needs an `RNA` assay with raw counts, the clustering column `RNA_snn_res.0.1` (clusters 0 to 3), and a `umap` reduction. | Built from the Sun et al. raw data: CNCB-NGDC GSA CRA001755 (BioProject PRJCA001543) |

The DE table used by Fig 4C-E and S5D-F (`DE_genes_sig.tsv`) is not an input. It is written by the Fig 4A/S5A script to `results/Fig4A_S5A_differential_expression/`.

## Notes

- The column order of `counts_matrix.tsv` matters. The scripts label columns 1 to 4 as MHC II+ and 5 to 8 as MHC II-.
- The subset annotation in Fig 4C-E is fixed by cluster ID in `SUBSET_MAP`. If you rebuild `mouse_mk.rds`, check `qc_cluster_marker_means.csv` after the first run and update the map if the cluster IDs have changed.
- To keep the inputs somewhere else, set `MK_BULK_COUNTS` and `MK_SC_RDS` instead of copying the files here.
