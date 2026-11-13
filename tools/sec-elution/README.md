# Installed SEC elution qualification

This fixed qualification exercises the public installed `measure` and
`measure.sec` APIs for SEC concentration units, detector-delay alignment,
serialized steps, no-op normalization, and a downstream LALS result. The
independent standard-library Python oracle is maintained separately; this R
producer records observed package outputs and does not calculate an expected
numerical oracle.

## Run

Install the candidate `measure` and `measure.sec` packages into the library
used by R, then run the producer with an output directory that does not
already exist:

```sh
Rscript --vanilla tools/sec-elution/qualify.R /tmp/sec-elution-qualification-run
```

The producer reads the frozen JSON inputs from its own folder. You can copy the
complete `tools/sec-elution` folder and invoke `qualify.R` by absolute path
from any working directory, including a path containing spaces:

```sh
Rscript --vanilla "/path with spaces/sec-elution/qualify.R" /tmp/sec-elution-qualification-run
```

The producer checks for `jsonlite` before reserving the output directory,
because it needs `jsonlite` to write failure receipts. If `jsonlite` is
missing, it reports a preflight failure and does not create the output root.
After reservation, failures such as missing `digest`, `recipes`, or installed
SEC packages leave `qualification-progress.json` and `run-status.json` with an
error. The runner also requires installed
`measure` and `measure.sec` namespaces with `Meta/package.rds`. It checks that
each loaded namespace path matches `find.package()` and records both paths and
package versions. A checkout HEAD describes only the producer checkout; it is
kept separate from installed package provenance and is null for a copied
standalone folder.

The frozen `fixtures.json` and `declaration.json` inputs are SHA-256 checked
before execution. The runner reserves the output root once, writes
`qualification-progress.json` and `run-status.json` while it works, and
retains both with `complete: false` and an error on failure. A completed run
contains `results.json`, `runtime.json`, and a status receipt with
`complete: true`, `stage: "complete"`, and `error: null`. Every run needs a new
output directory.

The mass controls are independent software cohorts: one recipe is trained on
the three declared training curves per case, then each current bake row is
retained from the run using that fixture's declared injected mass. This does
not model a per-row mass API or a single physical injection. The LALS values,
flow, delay, optical constant, and calibration constant are synthetic
declarations.

After production, run the independent standard-library verifier:

```sh
python3 tools/sec-elution/verify.py /tmp/sec-elution-qualification-run
```

For a copied folder, invoke its `verify.py` by absolute path with the same
output root. The verifier creates its independent receipt exclusively and
refuses to overwrite a receipt from an earlier run.

This qualification checks software numerical, unit, alignment, missingness,
and serialization behavior only. It does not establish detector calibration,
actual flow accuracy, sample recovery, baseline suitability, instrument
performance, multi-angle or radius-of-gyration capability, vendor interchange,
or generic operation-aware metadata preservation.
