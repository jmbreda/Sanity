""" Sort a MTX file by row and save it into *gzipped* sorted MTX file. This script use *sort* utility and expected to work in Lnux/OSX."""

import argparse
import gzip
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Sort a MTX file by row.")
    parser.add_argument("-i", "--input-file",
                        dest="input_file",
                        help="Input MTX file",
                        required=True)
    parser.add_argument("-o", "--output-file",
                        dest="output_file",
                        help="Output MTX file. Default name is 'sorted_' + input filename + '.gz'",
                        default = None)

    args = parser.parse_args()
    
    input_file = args.input_file
    if args.output_file is None:
        filename = os.path.basename(input_file)
        output_file = os.path.join(os.path.dirname(input_file), "sorted_" + filename + ".gz")
    else:
        output_file = args.output_file
    temp_file = output_file + "_tmp"

    # Read the input MTX file
    if input_file.lower().endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(input_file, "rt") as fin, \
        gzip.open(output_file, mode="wt", compresslevel=6) as fout:
        header = []
        start_sort = False
        for line in fin:
            if not start_sort:
                if line.startswith("%") or line.strip() == "":
                    header.append(line)
                elif len(line.split()) == 3:
                    header.append(line)
                    start_sort = True
                    # Write the header lines to the output file
                    fout.write("".join(header))
                    break
        
        with open(temp_file, mode="wt") as fout_temp:
            sort_process = subprocess.Popen(["sort", "-k1,1n"], stdin=subprocess.PIPE, stdout=fout_temp,
                                    text=True, encoding="utf-8")
            for line in fin:
                sort_process.stdin.write(line)
            sort_process.stdin.close()
            sort_process.wait()

            if sort_process.returncode != 0:
                os.remove(temp_file)
                os.remove(output_file)
                raise RuntimeError(f"'sort' failed with return code {sort_process.returncode}")
            
    
        with open(temp_file, mode="rt") as fin_temp:
            for line in fin_temp:
                fout.write(line)
        os.remove(temp_file)


if __name__ == "__main__":
    main()