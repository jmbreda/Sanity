"""Compare complete inference outputs with a reference executable."""
import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def compare(first, second):
    names = {p.name for p in first.iterdir()} - {"sanity_command.txt"}
    assert names == {p.name for p in second.iterdir()} - {"sanity_command.txt"}
    for name in names:
        if name.endswith(".npy"):
            a, b = np.load(first/name, allow_pickle=False), np.load(second/name, allow_pickle=False)
            assert a.shape == b.shape and a.dtype == b.dtype, name
            assert np.isfinite(a).all() and np.isfinite(b).all(), name
            np.testing.assert_array_equal(a, b, err_msg=name)
        else:
            assert (first/name).read_bytes() == (second/name).read_bytes(), name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--candidate", type=Path,
                        default=Path(__file__).resolve().parents[1]/"bin/Sanity")
    args = parser.parse_args()
    args.reference, args.candidate = args.reference.resolve(), args.candidate.resolve()
    counts = np.random.default_rng(42).poisson(2, (12, 17)).astype(float)
    counts[0] = 0
    counts[1] = 1
    counts[2] = np.geomspace(0.01, 1e6, 17)
    counts[3] = 0.05
    checks = 0
    with tempfile.TemporaryDirectory(prefix="sanity-inference-") as temporary:
        directory = Path(temporary)
        matrix = directory/"counts.tsv"
        with matrix.open("w") as stream:
            stream.write("\t".join(f"cell_{c}" for c in range(17))+"\n")
            for gene, row in enumerate(counts):
                stream.write(f"gene_{gene}\t"+"\t".join(f"{v:.17g}" for v in row)+"\n")
        sparse = directory/"counts.mtx"
        with sparse.open("w") as stream:
            stream.write(f"%%MatrixMarket matrix coordinate real general\n12 17 {np.count_nonzero(counts)}\n")
            for row, column in zip(*np.nonzero(counts)):
                stream.write(f"{row+1} {column+1} {counts[row,column]:.17g}\n")
        genes, cells = directory/"genes.txt", directory/"cells.txt"
        genes.write_text("".join(f"gene_{g}\n" for g in range(12)))
        cells.write_text("".join(f"cell_{c}\n" for c in range(17)))
        for source in [matrix, sparse]:
            for method in ["MAP", "MLE", "EAP", "MARG"]:
                for configuration in [[], ["-no_norm", "true", "-nbin", "19", "-vmin", "0.01", "-vmax", "20"]]:
                    formats = [[], ["-npy"]]
                    if method == "MAP":
                        formats += [["--bonsai-output-mode"]]
                    for output_format in formats:
                        pair = directory/f"pair_{checks}"
                        pair.mkdir()
                        for label, binary in [("before", args.reference), ("after", args.candidate)]:
                            command = [str(binary), "-f", str(source), "-d", str(pair/label),
                                       "-n", "2", "-e", "true", "-v_m", method,
                                       "-mtx_genes", str(genes), "-mtx_cells", str(cells),
                                       *configuration, *output_format]
                            process = subprocess.run(command, capture_output=True, text=True)
                            if process.returncode:
                                raise RuntimeError(process.stdout+process.stderr)
                        compare(pair/"before", pair/"after")
                        checks += 1
    print(f"PASSED: {checks} inference comparisons with exact arrays, text and identifier order")


if __name__ == "__main__":
    main()
