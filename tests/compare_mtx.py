"""Compare MTX inference with TSV input, including reader fallbacks.

Run after compiling: python3 compare_mtx.py
"""
import gzip
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def main():
    binary = Path(__file__).resolve().parents[1] / "bin" / "Sanity"
    generator = np.random.default_rng(901)
    counts = generator.poisson(generator.lognormal(0, 1, (48, 32)))
    counts[0, :] = 0
    with tempfile.TemporaryDirectory(prefix="sanity-mtx-") as temporary:
        directory = Path(temporary)
        genes = directory / "genes.txt"
        cells = directory / "cells.txt"
        genes.write_text("".join(f"gene_{index}\n" for index in range(48)))
        cells.write_text("".join(f"cell_{index}\n" for index in range(32)))

        def run(matrix, destination, threads=None):
            command = [str(binary), "-f", str(matrix), "-d", str(destination),
                       "-mtx_genes", str(genes), "-mtx_cells", str(cells),
                       "--bonsai-output-mode"]
            if threads is not None:
                command += ["-n", str(threads)]
            process = subprocess.run(command, capture_output=True, text=True)
            if process.returncode:
                raise RuntimeError(process.stdout + process.stderr)
            return process.stdout + process.stderr

        for fractional in [False, True]:
            values = counts.astype(float) / (4 if fractional else 1)
            table = directory / "counts.tsv"
            with table.open("w") as stream:
                stream.write("\t".join(f"cell_{index}" for index in range(32)) + "\n")
                for index, row in enumerate(values):
                    stream.write(f"gene_{index}\t" + "\t".join(map(str, row)) + "\n")
            reference = directory / f"reference_{fractional}"
            run(table, reference, 1)
            lines = ["%%MatrixMarket matrix coordinate real general",
                     f"48 32 {np.count_nonzero(values)}"]
            for row, column in zip(*np.nonzero(values)):
                lines.append(f"{row + 1} {column + 1} {values[row, column]:.17g}")
            contents = ("\n".join(lines) + "\n").encode()
            for compressed in [False, True]:
                matrix = directory / ("counts.mtx.gz" if compressed else "counts.mtx")
                matrix.write_bytes(gzip.compress(contents) if compressed else contents)
                for threads in [None, 1, 2, 12]:
                    output = directory / f"output_{fractional}_{compressed}_{threads}"
                    log = run(matrix, output, threads)
                    expected_reader = "existing" if fractional or compressed else "fast"
                    assert f"reader={expected_reader}" in log, log
                    if expected_reader == "fast":
                        assert f"up to {threads or 4} threads" in log, log
                    for field in ["delta", "d_delta", "mu", "variance"]:
                        np.testing.assert_array_equal(np.load(reference / (field + ".npy")),
                                                      np.load(output / (field + ".npy")))
                    for field in ["geneID.txt", "cellID.txt"]:
                        assert (reference / field).read_bytes() == (output / field).read_bytes()
    print("PASSED: 16 MTX runs match TSV inference exactly, including fractional/gzip fallbacks")


if __name__ == "__main__":
    main()
