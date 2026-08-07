
import glob
import numpy as np
import os
import shutil
import subprocess
import sys


def compare(dir1, dir2):
    files = glob.glob("{}/*".format(dir1))
    passed = 0

    for ctrlfile in files:
        x = os.path.split(ctrlfile)[1]
        if x == "sanity_command.txt":
            continue
        myfile = ctrlfile.replace(dir1, dir2)
        if not os.path.exists(myfile):
            print("File {} not found in {}".format(x, dir2))
            passed = 2
            return passed
        with open(ctrlfile, "rt") as fin1, \
            open(myfile, "rt") as fin2:
            files_printed = False
            for line1 in fin1:
                line1 = line1.strip("\n")
                line2 = fin2.readline().strip("\n")
                if line1 == line2:
                    continue
                data1 = line1.split("\t")
                data2 = line2.split("\t")
                for i in range(len(data1)):
                    if data1[i] != data2[i]:
                        with open("compare.log", "at") as fout:
                            if passed == 0:
                                passed = 1
                            if not files_printed:
                                fout.write("Files differ:\n{}\n{}\n".format(ctrlfile, myfile))
                                files_printed = True
                            numbers_close = np.isclose(float(data1[i]), float(data2[i]))
                            if not numbers_close:
                                passed = 2
                            fout.write("-----\n{}\n{}\tnumbers close: {}\n".format(data1[i],
                                                                                    data2[i],
                                                                                    numbers_close))
    return passed
                        

if __name__ == "__main__":
    sanity_out_dir = "sanity_output"

    methods = ["MAP", "EAP", "MLE", "MARG"]

    # Truncate the log once, before any method runs. compare() opens it in append mode as it finds
    # differences, so truncating inside the loop would discard the earlier methods' details.
    with open("compare.log", "wt"):
        pass

    failed_methods = []
    methods_with_differences = []
    for method in methods:

        default_results = "default_output_{}".format(method)

        print("Running Sanity comparison test for method {}".format(method))
        with open("compare.log", "at") as fout:
            subprocess.run(["../bin/Sanity", "-f", "count_table.tsv",
                "-d", sanity_out_dir,
                "-v_m", method,
                "-e", "1"], stderr=fout, stdout=fout)
        result = compare(default_results, sanity_out_dir)
        shutil.rmtree(sanity_out_dir)
        if result == 0:
            result = "PASSED"
        elif result == 1:
            # Values differ but are within numpy.isclose tolerances: treated as success.
            result = "PASSED with acceptable differences"
            methods_with_differences.append(method)
        else:
            result = "FAILED"
            failed_methods.append(method)
        print("Test for method {} passed: {}\n".format(method, result))

        if result != "PASSED":
            print("See compare.log for details")

    if failed_methods:
        print("FAILED for: {}. See compare.log for details.".format(", ".join(failed_methods)))
        sys.exit(1)
    if methods_with_differences:
        print("All methods passed, but the results were not exactly identical for: {}. "
              "The differences are within the numpy.isclose tolerances and count as passing; "
              "see compare.log for details.".format(", ".join(methods_with_differences)))
    else:
        print("All methods passed, with results exactly identical to the reference output.")
    sys.exit(0)
