# STOK — Knowledge-guided surrogate-assisted optimization for ecYeastGEM

This repository package is assembled from the authors' original **MTO-new** MATLAB research archive. It contains the STOK optimizer, the original ten condition-specific yeast enzyme-constrained model files, all ten stored knowledge bases, the original experiment scripts, data-derived variable mapping, and an available subset of figure-analysis resources.

**Release status: source/data organization and static verification only.** No MATLAB/COBRA/Gurobi execution has been performed in the packaging environment. This package is **not yet validated to reproduce the published figures**. See `docs/VALIDATION_REQUIRED.md` before citing it as a computationally reproduced release.

## Research question and objective

STOK tunes **4,261 model-effective catalytic-capacity parameters** in ecYeastGEM under 10 carbon/nitrogen conditions. It optimizes reference-growth deviation (`F1 = |mu_model(kcat)-mu_ref|`) and an ATP-associated reaction-flux criterion. These are **model-level objectives**, not measured intrinsic enzyme turnover constants, not direct experimental growth errors, and not total physiological ATP consumption.

## Files and provenance

- `src/stok/`: original STOK optimizer, surrogate-assisted search, similarity functions, and helpers.
- `src/Bio_test/`: original Bio experiment, source-model construction, and surrogate/knowledge-base procedures.
- `src/bio_objectives/bio1.m`—`bio10.m`: FBA objective evaluators. Bio1/Bio2 parameter-update loops have been restored with author's confirmation.
- `src/NSGAII/`, `src/NDSort/`, `src/Surr/`: helper algorithms and surrogate training used by source archive.
- `src/Database_source/`: generic benchmark offline-knowledge helpers.
- `models/`: original base ecGEM, ten condition-specific model files, 4,261 variable positions, C/N source files, and derived mapping CSV.
- `knowledge_base/benchmarks_bio/`: full ten originally uploaded archived knowledge `.mat` files (Git LFS tracked).
- `knowledge_base/benchmarks_generic/`: originally uploaded generic benchmark knowledge `.mat` files.
- `analysis/`: the uploaded analysis scripts and spreadsheets available in the work archive. Full historical `results/` (~2.3 GB) were NOT included, so reproducing every original figure requires those separate result files and any missing analysis code.
- `experiments/`: convenience entrypoint and preflight validation.
- `docs/`: provenance, settings, limitations and MATLAB reproduction checklist.

## Required MATLAB software

1. MATLAB (MToP original README advertises MATLAB >= R2020b); Statistics and Machine Learning Toolbox for `fitrgp` and Latin-hypercube sampling.
2. COBRA Toolbox and **Gurobi** MATLAB solver interface for `changeCobraSolver('gurobi','LP',0)` and `optimizeCbModel`.
3. **MToP v1.5**, separately obtained from its maintainers: https://github.com/intLyc/MTO-Platform . Its source and Git history are **not** republished here. `UniformPoint`, `W1_Optimum`, `W2_Optimum`, and `W3_Optimum` come from that platform (or corresponding MToP directories).

Exact MATLAB/COBRA/Gurobi versions and solver tolerances used for the original runs were not recoverable from the archive; authors should supply those before final scientific deposition.

## Run the included experiment script

Open MATLAB in the repository root and use:

```matlab
startup_STOK('C:\path\to\MTO-Platform');  % Windows example; edit location
verify_STOK_inputs();                  % checks required input files only
run_bio_experiments();                 % long and expensive: Bio1-Bio10 x 10 runs
```

`run_bio_experiments` wraps the original `build_pool_bio.m` script, which sets population size 50, FEsMax=1000 and 10 runs per task. The script uses `rng('shuffle')` per run, saves the initial RNG states, and writes result files under `results/Full/NSGAII/`. **This is a new execution of the packaged code, not guaranteed to regenerate the original historical runs**, which require their original RNG states and matched software versions.

**Important caveats before a full run:** `src/Bio_test/problem_family_bio.m` defines `lb`/`ub` using negative raw S-matrix coefficients. The main online `sas_ckt.m` currently does not directly use these bounds; other routines can. The source archive also contains a definition of `cooldown_FE=50` that is not enabled in the active phase-switch condition, and true function evaluation count can overshoot FEsMax due to batching. These issues are documented but **not silently changed** in the assembled project.

For details, see `docs/VALIDATION_REQUIRED.md` and `docs/REPRODUCIBILITY_SETTINGS.md`.

## Code/data reuse and third-party components

`src/*` was extracted from the submitted research archive. Authorship and redistribution permissions for individual helper files and model resources should be confirmed before public release. MToP, COBRA Toolbox and Gurobi are **external dependencies**, subject to their respective licenses. No blanket open-source license has been attached to this package.

## Publishing to GitHub

1. Review and test locally; start with a **private** repository.
2. Install Git and Git LFS (`git lfs install`). Several `.mat` and `.xlsx` files are tracked through `.gitattributes`.
3. Execute in the extracted root:

```bash
git init
git lfs install
git add .
git commit -m "Prepare STOK ecGEM reproducibility package"
git branch -M main
git remote add origin https://github.com/<YOUR_ACCOUNT>/STOK-ecGEM.git
git push -u origin main
```

4. Supply validation versions, solver options, run records and the complete figure-reproduction pipeline before publication. Create a fixed release and archive a citable DOI only after scientific validation.
