# Tests

This directory contains a simple test for Sanity.

## Overview

The test runs Sanity on a small count matrix (`count_table.tsv`) with the following options:
- `-v_max 0` or `1`
- `-e 1`

The test prints out results, which may vary:
- Results could be identical
- Results could differ but with numbers being close
- Results could differ significantly

## Requirements

Python script `compare.py` requires `numpy` library.

## Running the Test

1. **Compile Sanity** in the `src` project directory
2. **Navigate to the tests directory**:
    ```bash
    cd tests
    ```
3. **Run the test script**
    ```bash
    python compare.py
    ```
    Differences are saved in `compare.log` file.
