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
    Details on differences are saved in `compare.log` (overwritten on each run, not tracked in git).

Note: `compare.py` always exits with code 0, even when a method fails. Check the printed "PASSED"/"FAILED" lines (and `compare.log` for details) rather than the exit code.
