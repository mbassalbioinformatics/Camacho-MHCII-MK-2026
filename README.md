# Camacho-MHCII-MK-2026

Analysis code for:

> Camacho V, Wang KG, Hanč P, Carminita E, Becker IC, Lee DH, Bassal MA, Falchetti M, Maggi J, Barrachina M, von Andrian U, Carrascal M, Gautam D, Weng C, Sankaran VG, Italiano JE, Machlus KR. **MHC II-expressing bone marrow megakaryocytes are noncanonical antigen presenting cells and activate CD4+ T cells ex vivo.** Journal and citation to be added on publication.

This repository contains the code for the bulk RNA-seq analyses of sorted MHC II+ and MHC II- bone marrow megakaryocytes (MKs), including their integration with the Sun et al. (2021) single-cell MK atlas.

## Figure-to-code map

| Panel(s) | Analysis | Code |
|---|---|---|
| Fig 4A, Fig S5A | Differential expression, MHC II+ vs MHC II- MKs (DESeq2): volcano plot and top 50 heatmap | [`figures/Fig4A_S5A_differential_expression`](figures/Fig4A_S5A_differential_expression) |
| Fig S5B | Gene Ontology biological process enrichment (clusterProfiler) | [`figures/FigS5B_go_enrichment`](figures/FigS5B_go_enrichment) |
| Fig 4B, Fig S5C | Transcription factor activity (VIPER with DoRothEA) and TF networks | [`figures/Fig4B_S5C_tf_activity`](figures/Fig4B_S5C_tf_activity) |
| Fig 4C to 4E, Fig S5D to S5F | Integration with the Sun et al. (2021) single-cell MK atlas: unified panel heatmaps, correlation mapping, validation control, signature scoring | [`figures/Fig4C-E_S5D-F_sun_atlas_integration`](figures/Fig4C-E_S5D-F_sun_atlas_integration) |

Other analyses in the paper, including high-dimensional flow cytometry, proteomics and immunopeptidomics, were performed with the tools described in the STAR Methods and are not part of this repository.

## Layout

```
R/             helpers shared by the bulk RNA-seq scripts
figures/       one folder per figure panel group, each with its own README and script
data/          input data, not tracked in git (see data/README.md)
environment/   recorded R session information
results/       outputs written by the scripts, not tracked in git
run_all.R      installs any missing packages, then runs every analysis in order
setup.R        installs the required packages only (optional)
```

## Quick start

```bash
git clone https://github.com/mbassalbioinformatics/Camacho-MHCII-MK-2026.git
cd Camacho-MHCII-MK-2026
# place the input files in data/ (see data/README.md)
Rscript run_all.R
```

`run_all.R` first installs any missing R packages. It then runs the differential expression, because that writes the DE table used by the single-cell integration, followed by the other analyses. Each figure script can also be run on its own from the repository root, and writes to `results/<figure folder>/`.

## Software

Tested with R 4.6.0 and Bioconductor 3.23 on macOS (Apple silicon), using the CRAN and Bioconductor packages listed in [`R/required_packages.R`](R/required_packages.R). `run_all.R` installs any that are missing. To install them without running the analyses:

```bash
Rscript setup.R
```

Every script records its own `sessionInfo.txt` in its results folder. The tested environment for each analysis is in [`environment/`](environment).

## Data availability

- **Bulk RNA-seq** of sorted MHC II+ and MHC II- MKs: NCBI GEO, accession **GSE328415**.
- **Immunopeptidomics:** ProteomeXchange, accession PXD076501.
- **Reference single-cell atlas:** Sun S, et al. Single-cell analysis of ploidy and the transcriptome reveals functional and spatial divergency in murine megakaryopoiesis. *Blood* 2021;138(14):1211-1224. [doi:10.1182/blood.2021010697](https://doi.org/10.1182/blood.2021010697). Raw data: CNCB-NGDC BioProject [PRJCA001543](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001543) (GSA CRA001755).

## Citation

If you use this code, please cite the paper above and the archived release of this repository (see [`CITATION.cff`](CITATION.cff)). The Zenodo DOI will be added here on release.

## License

MIT. See [`LICENSE`](LICENSE).

## Questions

Please open an issue on this repository.
