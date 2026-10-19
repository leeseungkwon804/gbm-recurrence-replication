# Recurrent glioblastoma expression predictable from the primary tumour

Analysis code for:

> **Recurrent Glioblastoma Expression Is Predictable from the Primary Tumour in Two
> Independent Cohorts.** *Diagnostics* (submitted 2026).

## What this repository contains

All code used to produce every number, table and figure in the paper. No expression
value is combined across cohorts anywhere in the pipeline: each patient's
`recurrent − primary` difference is computed inside its own dataset, and only those
differences are compared.

```
scripts/      analysis pipeline (a1–a3 cohort preparation, gb1–gb2 gene-level models,
              vd1 validation, ep/ external-platform projection, tn/ time-to-relapse
              and permutation nulls, fig_analyses.py figure inputs)
figs/         figure scripts (fig1–fig4, figS, shared style.py)
```

## Data

The raw data are **not redistributed here.** Both cohorts are public:

| Cohort | Accession | Platform |
|---|---|---|
| A — Taipei Medical University, Taiwan | [GSE339484](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339484) | bulk RNA-seq, Illumina NovaSeq 6000 (GPL24676) |
| B — University of Basel, Switzerland | [GSE254874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254874) | NanoString nCounter PanCancer IO 360 |

Download both series from GEO into `data/` before running. Gene symbol mapping uses
[mygene.info](https://mygene.info) (queried 15 September 2026).

## Reproducing

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy pandas statsmodels scikit-learn matplotlib mygene
python3 scripts/a1.py          # cohort A preparation
python3 scripts/a2.py          # cohort B preparation
python3 scripts/a3.py          # within-patient differences
python3 scripts/gb1.py         # gene-level models
python3 scripts/vd1.py         # cross-cohort replication
python3 figs/fig1.py           # ... fig2–fig4, figS
```

Seeds are fixed in each script; permutation nulls use the seed recorded in the paper.

## Licence

MIT (see `LICENSE`). Please cite the paper if you use this code.
