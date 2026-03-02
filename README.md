# sc-HCC-AFP

Code repository for a hepatocellular carcinoma study integrating single-cell RNA-seq, spatial transcriptomics, bulk RNA-seq, and clinical data to investigate CXCL10+ DC-Treg interactions in AFP-positive tumors.

## Manuscript

Manuscript in preparation. Publication and citation details will be added after formal release.

## Figure Scripts

### Main figures
- Figure 1: `Fig1.R`, `Fig1_d-e_HCC_scanpy.py`
- Figure 2: `Fig2.R`
- Figure 3: `Fig3.R`
- Figure 4: `Fig4.R`
- Figure 5: `Fig5.R`
- Figure 6: `Fig6.R`

### Supplementary figures
- Figure S1: `FigS1.R`
- Figure S2: `FigS2.R`
- Figure S3: `FigS3.R`
- Figure S4: `FigS4.R`
- Figure S5: `FigS5.R`

### Spatial mapping
- CellTrek processing: `CellTrek.R`

## Minimal Setup

- R (recommended: >= 4.2)
- Python (recommended: >= 3.9)
- Core R packages used across scripts include: `Seurat`, `harmony`, `GSVA`, `ComplexHeatmap`, `clusterProfiler`, `survival`, `survminer`
- Core Python packages used in Scanpy analysis include: `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`

## Usage Notes

- Place required inputs under `data/` using the filenames expected in each script.
- Run scripts per figure as needed.
- Some scripts contain fixed local paths (`setwd(...)`); update them to your environment before execution.

## License

MIT License. See `LICENSE`.
