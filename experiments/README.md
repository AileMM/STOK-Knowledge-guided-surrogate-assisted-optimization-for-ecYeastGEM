# Experiment entrypoints

`run_bio_experiments.m` wraps the uploaded Bio driver. From the root in MATLAB:

```matlab
startup_STOK('C:\path\to\MTO-Platform');
verify_STOK_inputs();
run_bio_experiments();
```

The original script chooses `NSGAII` by `a=1`, nominal FEsMax=1000 and `runs=10`. It stores output under `results/Full/NSGAII/` (gitignored). Before scientific claims of full reproducibility, verify exact FE consumption and missing comparator configurations as described in `docs/VALIDATION_REQUIRED.md`.

The original `src/Bio_test/generate_bio_SourceTasks.m` and `src/Database_source/*` are also included for provenance/algorithm examination. They were not transformed into guaranteed working stand-alone scripts.
