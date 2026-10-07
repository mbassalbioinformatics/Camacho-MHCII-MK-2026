# Reference output

The figures and tables produced by running `Rscript run_all.R` on the data described in [`data/README.md`](../data/README.md). They are provided so a re-run can be checked against them without regenerating anything.

Generated with R 4.6.0 and Bioconductor 3.23 on macOS (Apple silicon). The environment for each analysis is in [`environment/`](../environment).

## Contents

| Folder | Published panels | Files |
|---|---|---|
| [`Fig4A_S5A_differential_expression`](Fig4A_S5A_differential_expression) | Fig 4A, Fig S5A | `Fig4A_volcano.pdf`, `FigS5A_top50_DE_heatmap.pdf`, `DE_all_genes.tsv`, `DE_genes_sig.tsv` |
| [`FigS5B_go_enrichment`](FigS5B_go_enrichment) | Fig S5B | `FigS5B_GO_BP_dotplot.pdf`, `GO_BP_MHCpos.tsv`, `GO_BP_MHCneg.tsv` |
| [`Fig4B_S5C_tf_activity`](Fig4B_S5C_tf_activity) | Fig 4B, Fig S5C | `Fig4B_TF_activity_barplot.pdf`, `FigS5C_TF_networks_MHCpos.pdf`, `FigS5C_TF_networks_MHCneg.pdf`, `TF_activity_summary.tsv`, `TF_activity_per_sample.tsv` |
| [`Fig4C-E_S5D-F_sun_atlas_integration`](Fig4C-E_S5D-F_sun_atlas_integration) | Fig 4C to 4E, Fig S5D to S5F | Panel PDFs, `bulk_vs_subset_spearman.csv`, `bulk_nearest_subset.csv`, `mapping_validation_control.csv`, intermediate tables, annotation QC |

## Agreement with the published analysis

- **Differential expression (Fig 4A, S5A):** the same 753 genes (518 higher in MHC II+, 235 higher in MHC II-), with identical statistics to floating-point precision.
- **TF activity (Fig 4B, S5C):** identical delta activity for all 272 TFs. The top 25 and the network TFs match the published panels.
- **Single-cell integration (Fig 4C to 4E, S5D to S5F):** the correlation, nearest-subset and validation-control tables are identical to the published analysis, including 12 of 12 for the validation control.
- **GO enrichment (Fig S5B):** the input gene lists are identical. GO annotations are updated with each Bioconductor release, so the terms here are similar but not identical to the published panel, which was generated with an earlier annotation release.

## Checking a re-run

Statistics should match to floating-point precision for the same package versions. Figure layout can differ slightly between graphics devices and package versions.

```bash
Rscript run_all.R
diff results/Fig4A_S5A_differential_expression/DE_genes_sig.tsv \
     reference_output/Fig4A_S5A_differential_expression/DE_genes_sig.tsv
```
