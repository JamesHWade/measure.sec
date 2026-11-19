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
