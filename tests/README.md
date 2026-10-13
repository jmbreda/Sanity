# Tests

This directory contains a simple test for Sanity.

## Overview

The test runs Sanity once for each `-v_m` method (MAP, EAP, MLE, MARG), each time with `-e 1`, on a small count matrix (`count_table.tsv`). Each run's output is compared against the reference results in the matching `default_output_<METHOD>` folder.

The comparison for each method can come out as:
- Identical results
- Results that differ but are numerically close (checked with `numpy.isclose`, using its default tolerances)
- Results that differ significantly

## Requirements

Python script `compare.py` requires the `numpy` library (`pip3 install numpy` if you don't have it).

## Running the Test

1. **Compile Sanity** in the `src` project directory
2. **Navigate to the tests directory**:
    ```bash
    cd tests
    ```
3. **Run the test script**
    ```bash
    python3 compare.py
    ```
    Details on differences are saved in `compare.log`, which is truncated at the start of each run and covers all four methods. It is not tracked in git.

Note: `compare.py` exits with code 0 if every method passed, and 1 if any method FAILED, so it can be used in a script or CI. A method whose results differ but are within the `numpy.isclose` tolerances is reported as "PASSED with acceptable differences" and counts as a pass.

## Matrix Market reader

After compiling, run `python3 compare_mtx.py`. It compares complete Bonsai NPY
outputs and identifier order against TSV input using the default thread setting
and `-n 1`, `-n 2`, and `-n 12`. It covers the parallel integral-count path,
fractional counts, gzip input, and removal of an empty gene. Temporary inputs
and outputs are removed automatically.

Run `python3 compare_mtx_index.py` to check every index entry and total,
partition/buffer boundaries, the thread limit, reader fallbacks and malformed
input rejection. It compiles a small test driver using `CXX` (default `g++`),
requires NumPy and removes its temporary files automatically.

## Per-gene input rows

Run `python3 compare_rows.py` to check exact counts for integer/fractional MTX,
scientific notation, large values, duplicates, gzip, TSV, reverse reads and
buffer boundaries. It also checks concurrent reader instances, limited thread
teams, stale row indexes, malformed entries and truncated rows. The script
uses NumPy, compiles a temporary driver with `CXX` (default `g++`), and removes
its temporary files.

## Inference changes

Compile the previous version separately, then run:

```bash
python3 compare_inference.py --reference /path/to/previous/bin/Sanity
```

An optional `--candidate` selects a different candidate executable. This checks
all four variance methods on TSV and fractional MTX inputs, a broad count range,
empty genes, a custom variance grid and disabled cell-size normalisation.
Extended text, float64 arrays, Bonsai arrays and identifier order must match
exactly. Command metadata is excluded because the executable paths differ.
Inputs and outputs are generated in a temporary directory and removed.
