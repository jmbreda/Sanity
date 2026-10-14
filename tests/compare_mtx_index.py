"""Check complete index contents, partition boundaries and reader fallbacks."""
import gzip
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np


def expected_index(content):
    lines = content.splitlines(keepends=True)
    offset = 0
    rows = cells = None
    entries = []
    for line in lines:
        fields = line.strip().split()
        if not fields or fields[0].startswith(b"%"):
            offset += len(line)
            continue
        if rows is None:
            rows, cells, count = map(int, fields)
        else:
            row, cell = int(fields[0])-1, int(fields[1])-1
            entries.append((row, cell, float(fields[2]), offset))
        offset += len(line)
    assert len(entries) == count
    libraries, totals = np.zeros(cells), np.zeros(rows)
    lengths, offsets = np.zeros(rows, dtype=np.uint64), np.full(rows, -1, dtype=np.int64)
    for row, cell, value, offset in entries:
        if lengths[row] == 0:
            offsets[row] = offset
        lengths[row] += 1
        totals[row] += value
        libraries[cell] += value
    retained = lengths > 0
    return (rows, int(retained.sum()), cells), {
        "offsets": offsets[retained], "entries": lengths[retained],
        "row_indices": np.flatnonzero(retained).astype(np.int32),
        "gene_totals": totals[retained], "cell_totals": libraries}


def main():
    root = Path(__file__).resolve().parents[1]
    fixtures = {
        "integral": b"%%MatrixMarket matrix coordinate real general\n% comment\n9 7 11\n1 1 3\n1 2 0\n1 4 8\n3 2 1e0\n3 2 4\n3 7 7e0\n5 1 2.0\n5 5 +6\n8 2 9007199254740980\n8 3 10\n9 7 3",
        "crlf": b"%%MatrixMarket matrix coordinate integer general\r\n4 3 4\r\n\r\n1 1 1\r\n1 2 2\r\n3 1 4\r\n4 3 8\r\n",
        "empty": b"%%MatrixMarket matrix coordinate integer general\n2 3 0",
        "split_gene": ("%%MatrixMarket matrix coordinate real general\n3 4 1000\n" + "".join(f"2 {i%4+1} {i%7}\n" for i in range(1000))).encode(),
        "fractional": b"%%MatrixMarket matrix coordinate real general\n3 4 3\n1 1 .5\n1 2 2.25\n3 4 3.125\n",
        "large_total": b"%%MatrixMarket matrix coordinate real general\n3 4 2\n1 1 9007199254740992\n2 1 1\n",
        "gzip": b"%%MatrixMarket matrix coordinate integer general\n3 4 2\n1 1 5\n3 4 6\n",
        "buffer_boundary": b"%%MatrixMarket matrix coordinate real general\n3 4 600003\n" + b"2 1 1234567890\n"*600003,
    }
    checks = 0
    with tempfile.TemporaryDirectory(prefix="sanity-index-") as temporary:
        directory = Path(temporary)
        driver = directory/"mtx_index"
        command = [*shlex.split(os.environ.get("CXX", "g++")), "-O2", "-std=c++17", "-fopenmp",
                   "-I"+str(root/"src"), str(root/"tests/mtx_index.cpp"),
                   str(root/"src/helpers/io/ReadInputFiles.cpp"), "-lz", "-o", str(driver)]
        subprocess.run(command, check=True)
        for name, content in fixtures.items():
            compressed = name == "gzip"
            matrix = directory/(name + (".mtx.gz" if compressed else ".mtx"))
            matrix.write_bytes(gzip.compress(content) if compressed else content)
            dimensions, expected = expected_index(content)
            for threads in [1, 2, 4, 12]:
                output = directory/f"{name}_{threads}"
                process = subprocess.run([str(driver), str(matrix), str(threads), str(output)],
                                         capture_output=True, text=True, check=True)
                fallback = name in ["fractional", "large_total", "gzip"]
                assert ("reader=existing" in process.stdout) == fallback, process.stdout
                assert tuple(map(int, (output/"dimensions.txt").read_text().split())) == dimensions
                for field, values in expected.items():
                    np.testing.assert_array_equal(np.fromfile(output/(field+".bin"), dtype=values.dtype), values)
                checks += 1
        limited = directory/"limited"
        environment = dict(os.environ, OMP_THREAD_LIMIT="2", OMP_DYNAMIC="FALSE")
        process = subprocess.run([str(driver), str(directory/"integral.mtx"), "12", str(limited)],
                                 env=environment, capture_output=True, text=True, check=True)
        assert "threads=2" in process.stdout
        for field, values in expected_index(fixtures["integral"])[1].items():
            np.testing.assert_array_equal(np.fromfile(limited/(field+".bin"), dtype=values.dtype), values)
        checks += 1
        malformed = {
            "unsorted": "2 1 2\n1 2 3\n", "column": "1 9 2\n2 1 3\n",
            "row": "4 1 2\n2 1 3\n", "negative": "1 1 -2\n2 1 3\n",
            "fields": "1 1 2 4\n2 1 3\n", "truncated": "1 1 2\n",
            "nan": "1 1 nan\n2 1 3\n", "infinity": "1 1 inf\n2 1 3\n",
        }
        for name, body in malformed.items():
            matrix = directory/(name+".mtx")
            matrix.write_text("%%MatrixMarket matrix coordinate real general\n3 4 2\n"+body)
            for threads in [1, 12]:
                process = subprocess.run([str(driver), str(matrix), str(threads), str(directory/(name+str(threads)))],
                                         capture_output=True, text=True)
                assert process.returncode != 0, name
                checks += 1
    print(f"PASSED: {checks} MTX index, fallback, thread-limit and malformed-input checks")


if __name__ == "__main__":
    main()
