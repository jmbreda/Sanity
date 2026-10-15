"""Compare complete dense counts before inference, including row/buffer boundaries."""
import gzip
import json
import os
import shlex
import tempfile
from pathlib import Path
import subprocess

import numpy as np

root = Path(__file__).resolve().parents[1]
temporary = tempfile.TemporaryDirectory(prefix="sanity-rows-")
directory = Path(temporary.name)
driver = directory/"read_rows"
subprocess.run([*shlex.split(os.environ.get("CXX", "g++")), "-O2", "-std=c++17", "-fopenmp",
                "-I"+str(root/"src"), str(root/"tests/read_rows.cpp"),
                str(root/"src/helpers/io/DenseRowReader.cpp"), str(root/"src/helpers/io/ReadInputFiles.cpp"),
                "-lz", "-o", str(driver)], check=True)
generator = np.random.default_rng(731)
counts = generator.poisson(2, (48, 64)).astype(float)
counts[0] = 0
counts[1] = 1
counts[2, 0] = 5
fixtures = []
for name in ["integer", "fractional", "scientific", "fractional_scientific", "real_dynamic_range", "large", "crlf", "no_final_newline", "gzip", "duplicate", "buffer_boundary", "tsv", "tsv_gzip"]:
    values = counts.copy()
    if name in ["fractional", "fractional_scientific"]: values /= 10
    if name == "real_dynamic_range":
        values = np.exp(generator.uniform(-20,30,values.shape))
        values[0]=0
    if name == "large": values[2, 0] = float(2**53+3)
    compressed = name in ["gzip", "tsv_gzip"]
    tsv = name.startswith("tsv")
    path = directory / (name + (".tsv" if tsv else ".mtx") + (".gz" if compressed else ""))
    if tsv:
        lines = ["\t".join(f"cell_{index}" for index in range(64))]
        lines += [f"gene_{index}\t"+"\t".join(f"{item:.17g}" for item in row) for index, row in enumerate(values)]
    else:
        entries = []
        for row, column in zip(*np.nonzero(values)):
            if name == "duplicate" and row == 2 and column == 0:
                entries.append(f"{row+1} {column+1} 999")
            token = f"+{values[row,column]:.17e}" if name in ["scientific", "fractional_scientific"] else f"{values[row,column]:.17g}"
            suffix = " "*(1024*1024+37) if name == "buffer_boundary" and row == 2 and column == 0 else ""
            entries.append(f"{row+1} {column+1} {token}{suffix}")
        lines = ["%%MatrixMarket matrix coordinate real general", f"48 64 {len(entries)}", *entries]
    content = ("\n".join(lines) + ("" if name == "no_final_newline" else "\n")).encode()
    if name == "crlf": content = content.replace(b"\n", b"\r\n")
    path.write_bytes(gzip.compress(content) if compressed else content)
    fixtures.append((name, path, "tsv" if tsv else "mtx", values[1:]))

checks = []
for name, path, format, expected in fixtures:
    for mode, threads in [("row", 1), ("reverse", 1), ("parallel", 2), ("parallel", 12)]:
        target = directory / (f"{name}_{mode}_{threads}.bin")
        command = [str(driver), mode, str(path), format, str(target), str(threads)]
        with target.with_suffix(".log").open("w") as logfile:
            subprocess.run(command, stdout=logfile, stderr=subprocess.STDOUT, check=True)
        observed = np.fromfile(target, dtype=np.float64).reshape(expected.shape)
        np.testing.assert_array_equal(observed, expected[::-1] if mode == "reverse" else expected)
        checks.append({"fixture": name, "mode": mode, "threads": threads, "exact": True})
environment = dict(os.environ, OMP_THREAD_LIMIT="2", OMP_DYNAMIC="FALSE")
path = directory / "limited.bin"
subprocess.run([str(driver), "parallel", str(fixtures[0][1]), "mtx", str(path), "12"],
               env=environment, stdout=subprocess.DEVNULL, check=True)
np.testing.assert_array_equal(np.fromfile(path,dtype=np.float64).reshape(counts[1:].shape), counts[1:])
checks.append({"thread_limit": 2, "requested": 12, "exact": True})
with (directory / "stale.log").open("w") as logfile:
    process = subprocess.run([str(driver), "stale", str(fixtures[0][1]), "mtx",
                              str(directory/"stale.bin"), "12"],
                             stdout=logfile, stderr=subprocess.STDOUT)
assert process.returncode != 0
assert "does not match indexed gene" in (directory / "stale.log").read_text()
checks.append({"stale_index": "rejected"})
for name,text in {"column":"1 3 1\n", "row":"2 1 1\n", "fields":"1 1\n",
                  "negative":"1 1 -2\n", "nan":"1 1 nan\n", "infinity":"1 1 inf\n",
                  "extra":"1 1 1 2\n", "count":"1 1 nonsense\n"}.items():
    matrix=directory/("invalid_"+name+".mtx");matrix.write_text(text)
    for mode in ["raw", "raw_parallel"]:
        with matrix.with_suffix("."+mode+".log").open("w") as logfile:
            process=subprocess.run([str(driver),mode,str(matrix),"mtx",str(directory/"invalid.bin"),"12"],
                                   stdout=logfile,stderr=subprocess.STDOUT)
        assert process.returncode != 0, name+" accepted"
        checks.append({"invalid":name,"mode":mode,"rejected":True})
matrix=directory/"truncated.mtx";matrix.write_text("1 1 1\n")
process=subprocess.run([str(driver),"raw_truncated",str(matrix),"mtx",str(directory/"invalid.bin"),"12"],
                       stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
assert process.returncode != 0 and b"Unexpected end" in process.stderr
checks.append({"truncated_row":"rejected"})
(directory / "verification.json").write_text(json.dumps({"passed": True, "checks": checks},indent=2)+"\n")
print("COUNT VERIFICATION PASSED", len(checks), "checks", flush=True)

temporary.cleanup()
