# Calibration of pooled CRISPR screens — v1.0.1

This package accompanies **Model-conditional inference separates robust conditional fitness from biological evidence in tumor–immune CRISPR screens**. It links screen fidelity, empirical-null inference, held-out calibration and coherent spike recovery to guide/culture robustness and evidence classification. The statistical runner implements the frozen six-culture Lawson contrasts. General workflow concepts can be adapted to other assays; such adaptations require a new design and qualification plan.

## Status and scope

Public code and derived-result package prepared for GitHub release v1.0.1 at https://github.com/seefreewind/immune-escape-fitness-decomposition. Original code is MIT licensed; project-generated numerical result tables are CC BY 4.0 (see DATA_LICENSE.md). Upstream data and annotation terms remain with their providers. Raw source data, original author workbooks, guide sequences and MSigDB membership files are excluded. Annotation labels and study-derived module assignments are analysis outputs, not redistributed upstream membership collections. Reproduction scope is described below.

## Installation and recorded versions

Use Python 3.12 and R 4.4.3. Install Python packages with `python -m pip install -r requirements.txt`. For R use `Rscript install_R_dependencies.R` with R 4.4.x/Bioconductor 3.20; edgeR must be 4.4.0. Package builds and BLAS can affect numerical reproduction: inspect `R_session_info.txt` and compare derived results. The dependency installer does not overwrite an existing environment automatically; use an isolated environment. The statistical runner requires pandas/NumPy; the remaining Python packages support archived downstream modules.

## Input contract

Tab-separated UTF-8 counts: `sgRNA`, `Gene`, then exactly six configured samples, comparator cultures first. Preserve the entire guide library, including rows later excluded from inference, because raw library exposure is retained. Counts must be nonnegative integers, guides unique, genes nonempty. The first two columns also supply guide-to-gene annotation. Config records condition labels, sample order and replicate provenance. Biological blocks and paired designs are not implemented by the six-unpaired-culture runner. A separate baseline/timepoint input is optional conceptually; exact B16 baseline/temporal/author compatibility flags are supplied as frozen ancillary evidence in this release.

Lawson raw counts: GEO GSE149933. Extract the six B16 T17 samples MUS_103, MUS_104, MUS_105, MUS_112, MUS_113, MUS_114 from the original integer guide counts and map guides using the source annotation. Do not use CPM or prefiltered gene counts. Set the `counts` path in a copy of `example_config/b16.json`; relative paths resolve relative to that config. No download is performed implicitly.

## Small software demonstration

From the release root:

```sh
python run.py demo --output demo_output
python -m unittest discover -s tests -v
python run.py audit --output evidence_audit
```

The demonstration checks exact finite-reference arithmetic on fixed toy values. It is not a new scientific simulation benchmark, has no inferential status and cannot produce manuscript gene counts. `audit` summarizes the frozen evidence table and its modality labels without rerunning discovery. Output directories must not already exist.

## Full B16 statistical reproduction

```sh
python run.py full --config example_config/b16.json --output preflight_B16 --preflight-only
python run.py full --config example_config/b16.json --output full_B16
```

The full command validates inputs/environment, refits predictive fidelity, and proceeds only after passing all mandatory criteria. It then runs 1,000 training nulls (including the 200-screen pilot), 100 held-out nulls and 60 coherent-spike screens, with complete normalization/dispersion/QL refitting. It evaluates six culture omissions and conditional guide omissions, then reconstructs B16 labels using explicitly frozen baseline/temporal/author and downstream annotation/evidence inputs. A complete regeneration of these ancillary layers from raw external data requires the original resources, matching versions and staged source files described in `SOURCE_INPUTS.tsv`; the lightweight runner does not regenerate all seven biological layers from raw data in one command. No new enrichment is run. These limitations distinguish full primary statistical reproduction from full raw-data-to-manuscript reproduction.

Expected frozen B16 outputs: 18,710 genes; 3,999 empirical FDR≤0.10; 2,675 also effect-qualified; 1,492 operational HC; 1,354 residual; Tier A/B/C=0/179/1,313. Compare rerun scores and component flags with `data/B16_gene_inference.tsv`; investigate differences rather than overwriting reference outputs. Model qualification includes additional Lawson-specific benchmarks retained in archived source and manuscript Methods. A statistically passing custom run alone is not the complete historical CALIBRATED_SIGNAL decision.

The full command is computationally intensive (1,160 complete genome-wide model refits plus predictive/omission fits). Use logs under `OUTPUT/11_logs`; duration depends on hardware. This revision tested software and frozen-output audits, not another full 1,000-null run.

## Outputs and seeds

`INPUT_VALIDATION.json`, per-screen `NULL_FIDELITY.tsv`, `EMPIRICAL_INFERENCE.tsv`, `CALIBRATION_SUMMARY.tsv`, `SPIKE_SUMMARY.tsv`, `RELEASE_ROBUSTNESS.tsv`, `B16_RECONSTRUCTED_CLASSIFICATION.tsv`, and `STATUS.json`. Fidelity failures produce no empirical P/FDR. Failed later calibration leaves candidate probabilities provisional. Monte Carlo training uses the +1 denominator 1,001; BH uses every eligible gene. Inference base seed is 20261019 + 100000 × contrast index in [4T1, B16, MC38]; held-out offset=20000, spike offset=40000, spike selection offset=50000. Downstream modules use 20261020.

## Derived data and source code

`data/` contains gene scores/P/FDR/HC, fidelity/calibration summaries, module assignments, candidate evidence and sample metadata. Third-party raw/annotation files must be obtained from their official repositories under source terms. No raw guide sequence annotation is redistributed. `CODE_PROVENANCE.json` identifies byte-identical original scientific scripts; their file-oriented stages require the original relative workspace structure. The seven numbered directories define the release entry points and stage-specific contracts. Software tests exercise probability ties, BH correction, invalid count rejection, demo separation and frozen counts. Original scientific results and thresholds remain unchanged.

## Archive and citation

Zenodo preservation is enabled through the GitHub integration before the first release. Use the DOI issued for the archived v1.0.1 snapshot once available. CITATION.cff and .zenodo.json provide creator and version metadata. This release preserves the documented partial reproduction scope; it does not claim a fresh full 1,160-simulation or end-to-end external-data reproduction.

`08_statistical_audit/` includes the frozen finite-reference/BH and common-grid diagnostics, fixed-caller comparison, seeds and original staged audit sources. These preserve the inference boundaries accompanying the main results.
