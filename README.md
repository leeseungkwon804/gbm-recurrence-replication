# gbm-recurrence-replication

Reanalysis of two patient-matched primary/recurrent glioblastoma cohorts measured on different
platforms, testing whether the expression change between the first and the second resection is shared
and reproducible across cohorts.

Manuscript: *Diagnostics*, manuscript ID diagnostics-4669893. Licence: MIT.

No new data were generated. Every input is public.

---

## 1. What you need before running anything

### 1.1 Software

Python 3.11.10 with the versions used for the published results:

| Package | Version |
|---|---|
| numpy | 2.4.6 |
| pandas | 2.3.3 |
| scipy | 1.17.1 |
| statsmodels | 0.15.0 |
| scikit-learn | 1.9.1 |
| matplotlib | 3.11.2 |

```
pip install -r requirements.txt
```

The figures use Arial. If Arial is not installed, matplotlib falls back to another sans-serif face and the
figures will render with different metrics; the numbers are unaffected.

### 1.2 Input files

Four inputs are required. Three are downloaded from public archives, one is included in the repository.

| File | Where it goes | Source |
|---|---|---|
| `GSE339484_counts_All_Samples.txt.gz` | `scripts/` | GEO series GSE339484, supplementary file (cohort A, bulk RNA-seq, 14 pairs) |
| `GSE254874_processed.txt.gz` | `scripts/` | GEO series GSE254874, processed matrix (cohort B, NanoString IO 360, 17 pairs) |
| `GSE339484_family.soft.gz` | `data/external/geo_meta/` | GEO series GSE339484, SOFT family file (per-sample grade and IDH annotation) |
| `h.all.v2024.1.Hs.symbols.gmt` | `data/external/` | MSigDB Hallmark collection v2024.1 |

Also in the repository and used directly:

| File | Contents |
|---|---|
| `scripts/GSE254874_relapse.tsv` | patient identifier and relapse interval in days, transcribed from the GEO sample records |
| `scripts/out/ens2sym.json` | Ensembl gene identifier to gene symbol mapping, queried from mygene.info on 15 September 2026 and frozen so that the mapping does not drift |
| `scripts/out/replication.tsv` | the 2026-09-15 per-gene result, used only as an assertion target so that the pipeline fails loudly if a later change alters an earlier number |

The mapping file is frozen deliberately. Re-querying mygene.info on 14 October 2026 changed three symbols,
none of them panel genes, and left every result unchanged; re-querying it today may change others.

---

## 2. Reproducing the principal analyses

```
python scripts/fig_analyses.py
```

This is the single entry point. It reads the four inputs above, re-executes the original pipeline, asserts
the previously recorded numbers before computing anything new, and writes every intermediate table used by
the manuscript to `data/derived/`. Run time is about 20 minutes on one core, dominated by the permutation
and bootstrap loops. Random seeds are fixed; hyper-parameters are fixed in advance and not tuned.

It asserts, and will stop if any of these fail:

* cohort A reduces to 23,084 expressed genes and cohort B to 784 panel genes, with 699 measured in both;
* 14 and 17 complete patient pairs are found;
* 67 genes reach q < 0.05 in cohort B across the 699 shared genes, 64 of them with the same sign in
  cohort A, and the 699-gene Spearman correlation rounds to 0.563;
* the per-gene table matches `scripts/out/replication.tsv` to floating-point tolerance.

### 2.1 Outputs and what uses them

| Output in `data/derived/` | Contents | Used by |
|---|---|---|
| `summary.json` | every headline number quoted in the manuscript | all figures, Table 2 |
| `TableS1_pergene_699.tsv` | per-gene results for the 699 shared genes, both cohorts | Figures 2, 3; Table S3 |
| `A_genomewide.tsv` | per-gene results for all 23,084 cohort-A genes | Figure 2a, Figure 4d |
| `B_panel784.tsv` | per-gene results for all 784 cohort-B panel genes | Figure 2b |
| `pca_A.tsv`, `pca_B.tsv`, `qc_A.tsv`, `qc_B.tsv`, `samplecorr_*.tsv` | sample structure and quality | Figure 1b–d, Figures S1, S2 |
| `cat.tsv` | concordance-at-the-top curves | Figure 3d |
| `permutation_nulls.tsv` | sign-flip and gene-label permutation nulls | Figures 3f, 3g, S6 |
| `bootstrap.tsv` | patient bootstrap draws | Figure S8 |
| `ml_roc.tsv`, `ml_lopo_predictions.tsv` | classifier curves and per-patient predictions | Figures 3e, S7 |
| `modules.json`, `module_sample_scores.tsv`, `TableS2_modules.tsv`, `S9_module_sensitivity.tsv` | module definitions, scores and sensitivity | Figure 4a–c, 4e–f; Table 2; Figure S10 |
| `gsea_hallmark.tsv` | Hallmark enrichment in both cohorts | Figure 4d, Figure S9 |
| `S5_normalisation.tsv` | replication under alternative normalisations | Figure S5 |
| `S7_pair_influence.tsv` | leave-one-pair-out influence | Figure S7 |
| `S8_filtered_panel_genes.tsv`, `S8_resampling_diag.json` | panel genes removed by the expression filter; resampling diagnostics | Figure S8 |
| `TableS3_pairs.tsv` | per-pair metadata | Table S1 |

## 3. Reproducing the figures

Run after `fig_analyses.py`. Each script writes a PDF and a 600 dpi PNG next to itself and prints the
panel-letter coordinates so that row and column alignment can be checked numerically.

```
cd figs
python fig1.py      # Figure 1  design, cohorts, sample structure, gene space
python fig2.py      # Figure 2  within-cohort change
python fig3.py      # Figure 3  replication, classifier, permutation nulls
python fig4.py      # Figure 4  modules and pathways
python figS.py      # Figures S1-S10
```

`figs/style.py` holds the shared style and the 12-column grid: one panel letter is one axes slot, every
row has a single height, and panel letters are placed in figure coordinates so that letters in a row share
a y coordinate and letters in a column share an x coordinate. Each script calls `check()`, which prints the
letter coordinates and reports any misalignment.

## 4. Mapping from manuscript item to code

| Manuscript item | Produced by | Reads |
|---|---|---|
| Table 1 (cohort characteristics) | hand-assembled from the GEO records and the original reports | `data/derived/TableS3_pairs.tsv` |
| Table 2 (module changes) | `scripts/fig_analyses.py`, module section | `summary.json` |
| Figure 1 | `figs/fig1.py` | `pca_*.tsv`, `qc_*.tsv`, `samplecorr_*.tsv`, `summary.json` |
| Figure 2 | `figs/fig2.py` | `A_genomewide.tsv`, `B_panel784.tsv`, `TableS1_pergene_699.tsv` |
| Figure 3 | `figs/fig3.py` | `TableS1_pergene_699.tsv`, `cat.tsv`, `ml_roc.tsv`, `permutation_nulls.tsv`, `summary.json` |
| Figure 4 | `figs/fig4.py` | `modules.json`, `module_sample_scores.tsv`, `gsea_hallmark.tsv` |
| Figures S1–S10 | `figs/figS.py` | the remaining `data/derived` tables listed in Section 2.1 |
| Section 3.3 replication statistics | `scripts/fig_analyses.py` | — |
| Section 3.5 classifier | `scripts/fig_analyses.py`, classifier section | — |
| Section 3.6 modules and pathways | `scripts/fig_analyses.py`, module and GSEA sections | — |

## 5. Analyses added during revision

These were added in response to peer review and are reported in the revised manuscript and its
supplementary tables. Each is a standalone script that reads only files already produced by
`fig_analyses.py` or the deposited inputs, and writes to `revision/`.

| Script | Question | Principal outputs |
|---|---|---|
| `revision/c1_locked_classifier.py` | how much of the cross-cohort classifier performance depends on standardising with the target cohort | `C1_classifier_settings.tsv`, `C1_roc_curves.tsv` |
| `revision/c2_15pair.py` | cohort-B eligibility, per-sample quality, and the principal analyses in every possible 15-pair subset | `C2_cohortB_sample_qc.tsv`, `C2_15pair_subsets.tsv`, `C2_cohortA_subsets.tsv` |
| `revision/c3_patient_uncertainty.py` | patient-level uncertainty, with the gene selection held fixed and repeated | `C3_uncertainty.json`, `C3_leave_one_pair_out.tsv` |
| `revision/c4_pathway_null.py` | pathway enrichment under a null that preserves gene correlation, and under a common gene universe | `C4_hallmark_full.tsv`, `C4_hallmark_common_universe.tsv` |
| `revision/c5_change_score.py` | projection and magnitude of the change vector against the relapse interval | `C5_change_scores.tsv`, `C5_per_patient_scores.tsv` |
| `revision/c6_fdr.py` | observed minimum adjusted p-value in each testing family | `C6_fdr_families.json` |
| `revision/c7_cohortA_preprocessing.py` | nature of the deposited cohort-A values and the origin of the column-sum difference | `C7_column_sum_summaries.tsv` |
| `revision/vecwilcox.py` | vectorised signed-rank test used by `c4_pathway_null.py` | module, not a script |

`revision/vecwilcox.py` reproduces `scipy.stats.wilcoxon` under its default settings, including the
handling of zero differences, the switch from the exact distribution to the normal approximation when
zeros or ties are present, and the absence of a continuity correction. `c4_pathway_null.py` re-checks this
against scipy on the observed data every time it runs, over all 20,351 cohort-A genes and all 784 cohort-B
panel genes, and stops if any p-value differs by more than 1e-9.

Two of these analyses are Monte Carlo and their counts move slightly between runs with different seeds.
In particular, the number of Hallmark sets reaching q < 0.05 genome-wide in cohort A under gene-set
permutation is 17 or 18 depending on the seed; sets near the threshold move between runs, which is why the
revised manuscript reports the unadjusted p-values alongside the counts.

## 6. Known issue to fix before release

`scripts/fig_analyses.py` currently reads the GEO SOFT file for cohort A from an absolute path on the
machine where it was first run:

```python
soft = gzip.open('<absolute path on the original machine>/GSE339484_family.soft.gz', 'rt').read()
```

The same file is in the repository at `data/external/geo_meta/GSE339484_family.soft.gz`. The line must be
changed to read it from there before the repository can be run by anyone else:

```python
soft = gzip.open(os.path.join(P, 'data', 'external', 'geo_meta', 'GSE339484_family.soft.gz'), 'rt').read()
```

This is the only absolute path in the analysis or figure code; the rest of the pipeline resolves every
path relative to the repository root.

## 7. Data availability

Cohort A: GEO accession **GSE339484** (Taipei Medical University, Taiwan; bulk RNA-seq, 14 patient-matched
pairs). Cohort B: GEO accession **GSE254874** (University of Basel, Switzerland; NanoString nCounter
PanCancer IO 360, 17 patient-matched pairs). GSE254875 is the SuperSeries containing GSE254874 and the same
patients; it is not used as a further cohort.

Hallmark gene sets: MSigDB v2024.1. Gene symbol mapping: mygene.info, queried 15 September 2026.

## 8. Design rule

Expression values are never pooled or compared across cohorts. Each patient's recurrent-minus-primary
difference is computed within that patient's own cohort, and only those differences, statistics derived
from them, or models fitted within one cohort and applied to the other after standardisation within that
cohort, cross between cohorts. Any change to the code should preserve this; `fig_analyses.py` is written so
that the two cohorts' expression matrices are never concatenated.

## 9. Citation

Please cite the manuscript. If you use the module definitions or the shared 699-gene set, note that the
gene set was fixed before testing as the panel genes measured in both cohorts, and that the 67-gene set
discussed in the manuscript is selected by significance in cohort B and is not a validated signature: it
ranges from 19 to 95 genes when a single patient pair is left out.
