# What the SEC fixtures measure

From the repository root, run:

```sh
Rscript tools/fixture-qualification/standards.R . /tmp/sec-fixtures
```

The script replays the twelve polystyrene standards and compares all 129,612
rows with the distributed data. It captures the signal before drift, after drift,
and after noise. It also measures the maximum on each recorded trace before and
after the same linear edge subtraction used in the calibration tutorial.

Read `standards-qualification.pdf` to compare peak times. In
`standards-fit-summary.csv`, two neighboring standards elute out of order when
the raw maximum defines the peak; three do after edge subtraction. The fitted
cubic increases over 42 of 199 intervals within the latter calibration range.
This is a discrepancy in a fixed synthetic example, not an acceptance limit.
`standards-stage-summary.csv` gives the twelve individual results, including
the noiseless response remaining in the right baseline window.

The replay evaluates generator definitions without running the commands that
write package data. Outputs go to the directory supplied on the command line.

## Compare the unknowns on the supported interval

Install the checkout, then run
`Rscript tools/fixture-qualification/unknowns.R . /tmp/sec-fixtures`.
The script compares the installed data with the source files, then checks the
public recipe against an independent calculation on the same time grid.
It also compares both implementations with a two-component mixture whose mass
fractions give known Mn, Mw, and Mz.

In `unknown-moment-qualification.csv`, the fraction of positive area retained
is the trapezoidal integral of the baseline-corrected signal (negative sampled
values set to zero) from 14.6367 to 16.2083 minutes, divided by the corresponding
integral over the recorded 4–22 minutes. Only 19.3%–23.0% remains. It is a
fraction of observed detector response, not a fraction of the assigned molecular
weight. The nominal generator mapping is evaluated only on that same interval,
so the comparison does not require extrapolation.
