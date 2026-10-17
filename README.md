# STOK: Knowledge-Guided Surrogate-Assisted Optimization for ecYeastGEM

This repository provides the MATLAB implementation of **STOK**, a knowledge-guided, surrogate-assisted multi-objective optimization framework for condition-specific **effective** $k_{cat}$ parameterization of *Saccharomyces cerevisiae* enzyme-constrained genome-scale metabolic models (ecGEMs).

STOK uses previously optimized conditions as an offline knowledge base, selects relevant source information, transfers selected decision-variable values, and combines knowledge transfer with surrogate-assisted optimization. The optimized quantities are model-effective enzyme-capacity parameters, **not** independently validated intrinsic biochemical turnover constants.

## Repository contents

| Location | Description |
| --- | --- |
| `src/stok/` | Main STOK optimizer, knowledge-transfer procedures, similarity calculations and surrogate-assisted search |
| `src/Bio_test/` | Bio-condition experiment driver and offline knowledge-base construction (`Knowledge_bio.m`, `saea_source_bio.m`, `build_pool_bio.m`) |
| `src/bio_objectives/` | `bio1.m`–`bio10.m` ecGEM objective functions |
| `src/NSGAII/`, `src/NDSort/`, `src/Surr/` | Optimization and surrogate-model helpers |
| `models/` | Base ecYeastGEM, condition-specific models, and coefficient-index data |
| `models/metadata/` | Mapping of all 4,261 decision-variable positions and recorded condition-specific media-bound changes |
| `experiments/` | Convenience scripts for building the offline knowledge base, checking inputs and starting the main experiment |
| `knowledge_base/benchmarks_bio/` | **Locally generated** Bio1–Bio10 offline knowledge archives; large `.mat` outputs are intentionally not committed |
| `knowledge_base/benchmarks_generic/` | Selected general multi-objective benchmark knowledge files |
| `analysis/` and `docs/` | Available analysis materials and implementation/reproducibility notes |

The supplied base model identifies itself as `ecYeastGEM_batch_v8.3.4` and has 8,144 reactions. The parameterization used here has **4,261 decision-variable positions**. Bio1–Bio10 represent ten carbon/nitrogen-source conditions; their model files are under `models/bio_source_models/`.

## Software requirements

- MATLAB, with Statistics and Machine Learning Toolbox (e.g., `lhsdesign`, `fitrgp`).
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) for FBA evaluation.
- Gurobi, configured as the COBRA Toolbox LP solver.
- [MToP / MTO-Platform](https://github.com/intLyc/MTO-Platform) (the original research archive used MToP v1.5) for external optimization utilities.

Install and configure these dependencies separately. Their source code and licenses are not part of STOK. Exact historical MATLAB, COBRA, and Gurobi versions/solver tolerances should be recorded if exact numerical reproduction is required.

## How to run

Clone this repository (or download and extract it), open MATLAB, and set its current folder to the **repository root**. Replace the example MToP path with its actual local installation directory.

```bash
git clone https://github.com/AileMM/STOK-Knowledge-guided-surrogate-assisted-optimization-for-ecYeastGEM.git
```

```matlab
% 1. Add STOK and MToP functions to the MATLAB search path.
startup_STOK('C:\path\to\MTO-Platform');

% 2. Generate the Bio1–Bio10 offline knowledge base (required once).
generate_bio_knowledge_base();

% 3. Check that the required model, objective and knowledge files exist.
verify_STOK_inputs();

% 4. Run the main STOK optimization experiments.
run_bio_experiments();
```

### Step 1: Rebuild the offline knowledge base

The large Bio knowledge-base `.mat` files are **not uploaded to GitHub**, because they are outputs of an offline optimization procedure and can be regenerated from the supplied source code and model inputs.

The call chain is:

```text
experiments/generate_bio_knowledge_base.m
  └── src/Bio_test/Knowledge_bio.m
        └── src/Bio_test/saea_source_bio.m
              └── knowledge_base/benchmarks_bio/bio1-M2-d4261-k10.mat
                  ...
                  knowledge_base/benchmarks_bio/bio10-M2-d4261-k10.mat
```

`Knowledge_bio.m` invokes the **Bio-specific** source optimizer `saea_source_bio.m` (not the separate generic benchmark function `saea_source.m`). It saves archived source solutions, objective values and identified key-variable information. The source-stage code sets `FEsMax = 100` for its offline optimization. The `-k10` filename suffix comes from the original configuration and does **not**, by itself, indicate ten independent source-stage runs.

**Generate the knowledge base before running the main optimizer:** `build_pool_bio.m` loads all ten archives when it starts, then excludes the target condition when choosing sources for each optimization task. The generation step may take substantial time because it performs repeated ecGEM/FBA evaluations.

### Step 2: Run the STOK main experiments

`experiments/run_bio_experiments.m` invokes the original **`src/Bio_test/build_pool_bio.m`** driver. The current driver specifies an initial population of 50, a nominal maximum of 1,000 true evaluations per online run, and ten independent runs for each target Bio condition. It writes run results under `results/Full/NSGAII/`, and records the run RNG state in the saved outputs.

You can inspect or run `build_pool_bio.m` directly after initializing the MATLAB search path and generating the knowledge base; the convenience wrapper is provided to run it from the correct working directory.

## Reproducibility notes

- **Offline knowledge is stochastic.** Regenerated archives can differ from the historical source-task archives. Running the pipeline from source is possible, but it is not a guarantee of byte-for-byte reproduction of the manuscript's original numerical results.
- **The supplied objective functions evaluate model-level quantities.** The growth objective is the deviation from the original ecGEM growth reference, not directly from experimental growth measurements. The ATP-related objective is a solver-returned ATP-associated reaction-flux quantity rather than a direct physiological total-ATP-consumption measurement.
- **Model and software settings matter.** Model input files, enzyme-capacity constraints, random states and solver versions can affect optimization results. See `docs/REPRODUCIBILITY_SETTINGS.md` and `docs/VALIDATION_REQUIRED.md` for the code-derived settings and remaining execution-level checks.
- **Archived results and comparisons.** The full historical experiment result directory and all external comparator implementations are not necessarily included in this code-focused repository. Consequently, availability of this repository alone does not establish exact regeneration of every manuscript figure or table.

The repository was organized and statically reviewed; the complete MATLAB/COBRA/Gurobi workflow has **not** been independently executed as part of the packaging process.

## Third-party attribution and use

The project uses or references MToP, the COBRA Toolbox, Gurobi and the ecYeastGEM/GECKO model ecosystem. Please follow each project's licensing and citation requirements. Some optimization helper functions originate from the author's wider research archive; review redistribution permissions before applying a repository-wide software license.

## Paper and contact

This code accompanies the STOK manuscript on similarity-guided transfer optimization for effective catalytic-capacity parameterization in ecYeastGEM. A publication citation and versioned archival DOI can be added once available.
