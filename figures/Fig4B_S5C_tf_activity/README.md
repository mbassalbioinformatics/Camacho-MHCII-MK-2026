# Fig 4B and Fig S5C. Transcription factor activity

Script: [`Fig4B_S5C_tf_activity.R`](Fig4B_S5C_tf_activity.R)

## What this analysis does

Transcription factor (TF) activity is inferred for each sample and compared between MHC II+ and MHC II- MKs.

## Steps

1. The DESeq2 model is rebuilt with [`R/bulk_rnaseq_helpers.R`](../../R/bulk_rnaseq_helpers.R), identical to Fig 4A, and expression is variance-stabilized (`vst`, `blind = TRUE`).
2. Activity is inferred per sample with VIPER (`decoupleR::run_viper`) using mouse DoRothEA regulons of confidence A to C, with a minimum regulon size of 5.
3. Mean activity is calculated per TF and group. Delta activity is the MHC II+ mean minus the MHC II- mean.
4. **Fig 4B.** The 25 TFs with the largest absolute delta activity.
5. **Fig S5C.** Regulatory networks for the five TFs most increased in each group, each showing its ten DoRothEA targets with the lowest DE adjusted P, coloured by log2 fold change. Network layouts are seeded so they reproduce.

## Outputs

Written to `results/Fig4B_S5C_tf_activity/`.

| File | Content |
|---|---|
| `Fig4B_TF_activity_barplot.pdf` | Fig 4B |
| `FigS5C_TF_networks_MHCpos.pdf`, `FigS5C_TF_networks_MHCneg.pdf` | Fig S5C |
| `TF_activity_per_sample.tsv` | VIPER activity score for every TF and sample |
| `TF_activity_summary.tsv` | Mean activity per group and delta activity, ranked |
| `sessionInfo.txt` | R environment for the run |

Published result: the TFs most increased in MHC II+ MKs are Irf4, Rfx5, Irf9, Pparg and Irf1, and those most increased in MHC II- MKs are Smad4, Tbx21, Rel, Bach1 and Cebpb.

## Running

```bash
Rscript figures/Fig4B_S5C_tf_activity/Fig4B_S5C_tf_activity.R
```

Input: `data/counts_matrix.tsv` (see [`data/README.md`](../../data/README.md)).

## Requirements

Tested with R 4.6.0, DESeq2 1.52.0, decoupleR 2.17.0, viper 1.46.0, igraph 2.3.4 and ggraph 2.2.2 (Bioconductor 3.23), with DoRothEA regulons from the dorothea package. The full environment is in [`environment/sessionInfo_Fig4B_S5C.txt`](../../environment/sessionInfo_Fig4B_S5C.txt).

## Reproducibility

This script reproduces the published TF activity exactly: all 272 TFs, with identical delta activity to floating-point precision. The top 25 match published Fig 4B.
