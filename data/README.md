# Input data

Input files go in this folder. None of them are tracked in git. Inputs for each analysis are listed below. Further sections will be added as the code for other figures is added.

## Figure 4C

| File | Content | Source |
|---|---|---|
| `mouse_mk.rds` | Seurat object of the Sun et al. (2021) bone marrow MK scRNA-seq. It needs an `RNA` assay with raw counts, the clustering column `RNA_snn_res.0.1` (clusters 0 to 3), and a `umap` reduction. | Built from the Sun et al. raw data: CNCB-NGDC GSA CRA001755 (BioProject PRJCA001543) |
| `counts_matrix.tsv` | Bulk RSEM expected counts, genes by 8 samples, tab-separated. Row names have the form `ENSMUSG..._Symbol`. Columns 1 to 4 are MHC II+ and columns 5 to 8 are MHC II-. | RSEM gene-level output of the eight bulk libraries; accession to be added on deposition |
| `DE_genes_sig.tsv` | DESeq2 results for the significant genes, with at least the columns `gene` and `sig` (`Up_MHCpos` or `Up_MHCneg`). | DESeq2 on `counts_matrix.tsv` (absolute log2 fold change > 1, adjusted P < 0.05) |

### Notes

- The column order of `counts_matrix.tsv` matters. The script labels columns 1 to 4 as MHC II+ and 5 to 8 as MHC II-.
- The subset annotation is fixed by cluster ID in `SUBSET_MAP`. If you rebuild `mouse_mk.rds`, check `results/Fig4C/qc_cluster_marker_means.csv` after the first run and update the map if the cluster IDs have changed.
- To keep the inputs somewhere else, set `MK_SC_RDS`, `MK_BULK_COUNTS` and `MK_DE_TABLE` instead of copying the files here.
