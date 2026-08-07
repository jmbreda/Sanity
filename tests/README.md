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
