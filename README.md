# Clustering 2,700 PBMCs, and how stable the clusters are

A single-cell RNA-seq workflow in Python (Scanpy, Leiden) on the public 10x Genomics PBMC 3k dataset, with a
subsampling check on cluster stability and a written account of annotation attempts that failed.

**Report:** [`report/report.pdf`](report/report.pdf) (source: `report/report.qmd`; numbers come from `results.json`).

## Result in brief
- 2,700 droplets, 2,434 cells after adaptive (MAD) quality control; 8 Leiden clusters at resolution 0.5, annotated as 7 cell types.
- Clusters agree moderately with the Scanpy tutorial's Louvain clusters (ARI 0.61); per-type overlap is high for B, dendritic, CD14 monocyte and NK cells.
- Re-clustering 80% subsamples gives mean ARI 0.75 (SD 0.08): the broad structure is reproducible, boundaries are not.
- Two early annotation methods were wrong (a CD16 monocyte cluster labelled megakaryocyte; CD8 T cells labelled NK). They are documented in the report. The final rule table was written after inspecting the data, so annotation is not blind.

## What this adds beyond the standard workflow
- **The clusters are checked for stability, which tutorials skip.** Re-clustering 80% subsamples gives mean ARI 0.75.
- **The annotation is checked against a classifier that never saw this data.** CellTypist (`Immune_All_Low`) names
  all 8 clusters consistently with the hand-written rule table: classical monocytes on the CD14 cluster,
  non-classical on the CD16 cluster, cytotoxic T on CD8 T, helper T on both CD4 T clusters, and B, NK and dendritic
  directly. Its confidence ranges from 96% (NK) down to 31% (CD8 T) - lowest on the cluster an earlier attempt of
  mine got wrong, which is independent support for calling that boundary soft.
  The ARI between the two labellings is 0.50, but that mostly measures a resolution gap: CellTypist used 40 labels
  where the rule table uses 7.
- **It reproduces.** Re-running the whole pipeline from `environment.yml` in a fresh environment reproduced all 150
  values in `results.json` exactly.

## Reproduce
```bash
micromamba create -y -f environment.yml             # pins the versions the report quotes
micromamba run -n scrna-pbmc3k python src/analysis.py   # downloads the public data, writes figures/ and results.json
quarto render report/report.qmd                    # R with jsonlite and knitr
```

## Provenance
Written and run with the assistance of Claude Code (Anthropic). I am responsible for checking and understanding the results. No part derives from university coursework or assessed material.
`report/apa.csl` is the APA style from the [CSL styles repository](https://github.com/citation-style-language/styles) (CC BY-SA 3.0).
