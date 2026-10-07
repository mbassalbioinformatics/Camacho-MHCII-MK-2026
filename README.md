# Camacho-MHCII-MK-2026

Analysis code for Camacho et al., a study of MHC class II-expressing bone marrow megakaryocytes (MKs). Manuscript under revision.

## Status

This repository is being assembled figure by figure. The code for **Figure 4C** is complete. Code for the remaining computational panels will be added as it is finalized. Figure numbering follows the submitted manuscript and will be updated to match the published version.

## Figure-to-code map

| Figure | Panel(s) | Analysis | Code | Status |
|---|---|---|---|---|
| 1 | D | High-dimensional flow cytometry UMAP of bone marrow progenitors (Spectre) | pending | To be added |
| 3 | C, D; S3B, S3C | Proteomics of MKs, platelets, blood plasma and bone marrow fluid | pending | To be added |
| 4 | A; S4A | Bulk RNA-seq differential expression, MHC II+ vs MHC II- MKs (DESeq2) | pending | To be added |
| 4 | S4B | Gene Ontology biological process enrichment | pending | To be added |
| 4 | B; S4C | Transcription factor activity (decoupleR VIPER with DoRothEA) and TF networks | pending | To be added |
| 4 | C | Bulk RNA-seq integrated with the Sun et al. (2021) single-cell MK atlas | [`figures/Fig4C`](figures/Fig4C) | **Available** |
| 5 | B to E; S5A to S5D | MHC II immunopeptidomics (GibbsCluster, NetMHCIIpan, PANTHER) | pending | To be added |

Panels produced in FlowJo, Imaris, Fiji or similar software, and the standard statistical comparisons in the remaining panels, are described in the manuscript Methods rather than in this repository.

## Layout

```
figures/       one folder per figure, each with its own README and script(s)
data/          input data, not tracked in git (see data/README.md)
environment/   R session information for each analysis
results/       outputs written by the scripts, not tracked in git
```

## Quick start

```bash
git clone https://github.com/<your-username>/Camacho-MHCII-MK-2026.git
cd Camacho-MHCII-MK-2026
# place the input files in data/ (see data/README.md)
Rscript figures/Fig4C/Fig4C_bulk_vs_sun_atlas.R
```

Every script runs from the repository root and writes to `results/<figure>/`.

## Software

Tested with R 4.6.0 on macOS (Apple silicon). Package versions for each analysis are listed in its figure README, and the full environment is recorded in [`environment/`](environment).

## Data availability

- **Bulk RNA-seq** of sorted MHC II+ and MHC II- bone marrow MKs (n = 4 per group): accession to be added on deposition.
- **Reference single-cell atlas:** Sun S, et al. Single-cell analysis of ploidy and the transcriptome reveals functional and spatial divergency in murine megakaryopoiesis. *Blood* 2021;138(14):1211-1224. [doi:10.1182/blood.2021010697](https://doi.org/10.1182/blood.2021010697). Raw data: CNCB-NGDC BioProject [PRJCA001543](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001543) (GSA CRA001755).

## Citation

The paper citation and the Zenodo DOI for this code will be added here on release.

## License

MIT. See [`LICENSE`](LICENSE).

## Questions

Please open an issue on this repository.
