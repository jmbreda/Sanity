
import glob
import numpy as np
import os
import shutil
import sys


def compare(dir1, dir2):
    files = glob.glob("{}/*".format(dir1))
    passed = 0

    for ctrlfile in files:
        x = os.path.split(ctrlfile)[1]
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
    default_results = "default_output0"

    with open("compare.log", "at") as fout:
        pass  # Clear log file

    print("Running Sanity comparison test 1")
    os.system("../bin/Sanity -f count_table.tsv -d {} -max_v 0 -e 1".format(sanity_out_dir))
    result = compare(default_results, sanity_out_dir)
    shutil.rmtree(sanity_out_dir)
    if result == 0:
        result = "PASSED"
    elif result == 1:
        result = "PASSED with acceptable differences"
    else:
        result = "FAILED"
    print("Test 1 passed: {}\n\n\n".format(result))

    if result != "PASSED":
        print("See compare.log for details")

    print("Running Sanity comparison test 2")
    default_results = "default_output1"
    os.system("../bin/Sanity -f count_table.tsv -d {} -max_v 1 -e 1".format(sanity_out_dir))
    result = compare(default_results, sanity_out_dir)
    shutil.rmtree(sanity_out_dir)
    if result == 0:
        result = "PASSED"
    elif result == 1:
        result = "PASSED with acceptable differences"
    else:
        result = "FAILED"
    print("Test 2 passed: {}".format(result))

    if result != "PASSED":
        print("See compare.log for details")