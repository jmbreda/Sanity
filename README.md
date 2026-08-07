# Sanity

**Sa**mpling **N**oise based **I**nference of **T**ranscription Activit**Y** : Filtering of Poisson noise on a single-cell RNA-seq UMI count matrix

Single-cell RNA sequencing normalization algorithm presented in the publication [Bayesian inference of gene expression states from single-cell RNA-seq data -
J Breda, M Zavolan, E van Nimwegen - Nature Biotechnology, 2021](https://www.nature.com/articles/s41587-021-00875-x).

Sanity infers the log expression levels *x<sub>gc</sub>* of gene *g* in cell *c* by filtering out 
the Poisson noise on the UMI count matrix *n<sub>gc</sub>* of gene *g* in cell *c*.



### Reproducibility
The raw UMI count and normalized datasets mentioned in benchmarking in the associated [publication](https://www.nature.com/articles/s41587-021-00875-x) are available on [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4009187.svg)](https://zenodo.org/record/4009187). Files are named `[dataset name]_UMI_counts.txt.gz` and `[dataset name]_[tool name]_normalization.txt.gz`.

The scripts used for running the benchmarked normalization methods and for making the figures of the preprint are in the reproducibility folder.

## Input

* `-f`: UMI count matrix: *(N<sub>g</sub> x N<sub>c</sub>)* matrix with *N<sub>g</sub>* the number of genes and *N<sub>c</sub>* the number of cells. Format: tab-separated, comma-separated, or space-separated values. (`'path/to/text_file'`)

| GeneID | Cell 1 | Cell 2 | Cell 3 | ... |
|:-------|:------:|:------:|:------:|:---:|
| Gene 1 | 1.0 | 2.0 | 0.0 | ... |
| Gene 2 | 6.0 | 3.0 | 1.0 | ... |
| ... | ... | ... | ... | ... |

The count matrix can be compressed with gzip (with `.gz` extension).

* `-f`: (Alternatively) Matrix Market File Format: Sparse matrix of UMI counts. Automatically recognized by `.mtx` extension of the input file. Example: `matrix.mtx` by cellranger 2.1.0 and 3.1.0 (10x Genomics). (`'path/to/text_file.mtx'`)
**Important**: MTX file should be sorted by row (gene) indices! If not then you can use `sort_mtx_by_row.py` script in the `Sanity/scripts` folder to sort it. See [Sorting MTX file by row (gene) indices](#sorting-mtx-file-by-row-gene-indices).
The MTX file can be compressed with gzip (with `.gz` extension).
	* `-mtx_genes`: (optional) Gene ID file: text file with one gene ID per line. The order of gene IDs should match the order of genes in the count matrix. Examples: `genes.tsv` by cellranger 2.1.0 and `features.tsv` by cellranger 3.1.0 (10x Genomics). (`'path/to/text_file'`)
	* `-mtx_cells`: (optional) Cell ID file: text file with one cell ID per line. The order of cell IDs should match the order of cells in the count matrix. Examples: `barcodes.tsv` by cellranger 2.1.0 and 3.1.0 (10x Genomics).  (`'path/to/text_file'`)
* `-d`: (optional) Destination folder (`'path/to/output/folder'`, default: `cwd`)
* `-n`: (optional) Number of threads (integer, default: `4`)
* `-e`: (optional) Print extended output (Boolean, `'true', 'false', '1'` or `'0'`, default: `false`)
* `-v_m`: (optional, expert-user-only) Choose the method to estimate gene-variances $v_g$. In the MAP, EAP, MLE-options, one value for $v_g$ is fixed, and the corresponding gene expression estimates are returned, in the MARG-option, the gene expression estimates are obtained by marginalizing $v_g$. The options are:
    * MAP (**default**): Use the maximum a posteriori estimate for $v_g$.
    * EAP: expected value of $v_g$ over the posterior.
    * MLE: maximum likelihood estimate for $v_g$.
    * MARG: (original method) The reported $v_g$ is the EAP estimate, but unlike the EAP option above, this value is not fixed when computing the gene expression estimates — instead, the gene expression estimates are obtained by marginalizing over the posterior of $v_g$.
* `-vmin/-vmax`: (optional, expert-user-only) Minimal and maximal considered values of the variance in log transcription quotients (double, default: *v<sub>min</sub>=*`0.001` *v<sub>max</sub>=*`50`)
* `-nbin`: (optional, expert-user-only) Number of bins for the variance in log transcription quotients (integer, default: `160`)
* `-no_norm`: (optional, expert-user-only) Option to skip cell size normalization (Boolean, `'true', 'false', '1'` or `'0'`, default: `false`)

## Output

* `log_transcription_quotients.txt`: This file contains the estimated values of the log-transcription quotients (LTQs) for each gene in each cell. The LTQ *x<sub>gc</sub>* of gene *g* in cell *c* corresponds to the estimated logarithm of the fraction of mRNAs in cell *c* that belong to gene *g*. The LTQs are thus normalized such that *&Sigma;<sub>g</sub> exp(x<sub>gc</sub>) = 1* for each cell *c*. In order to get an estimate of the number of mRNAs for gene *g* in cell *c* one would thus need to multiply *exp(x<sub>gc</sub>)* by the estimated total number of mRNAs *M* in the cell.

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ... |
  |:-------|:------:|:------:|:------:|:---:|
  | Gene 1 | -13.7227 | -13.722 | -13.729 | ... |
  | Gene 2 |  -9.96744 | -10.2522 | -10.1453 | ... |
  | ... | ... | ... | ... | ... |
  
* `ltq_error_bars.txt`: Table with the error-bars on the estimates of the LTQs *x<sub>gc</sub>* for each gene *g* in each cell *c*.

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ... |
  |:-------|:------:|:------:|:------:|:---:|
  | Gene 1 | 0.630111 | 0.630198 | 0.624802 | ... |
  | Gene 2 | 0.315551 | 0.325912 | 0.301861 | ... |
  | ... | ... | ... | ... | ... |


## Extended output (optional)

* `mu.txt`: Estimated average LTQ *&mu;<sub>g</sub>* of each gene *g* averaged over all cells.One value per line. Order corresponds to the order of genes in the `geneID.txt` file.
* `d_mu.txt`: Error bars on the inferred mean LTQs *&mu;<sub>g</sub>*. One value per line. Order corresponds to the order of genes in the `geneID.txt` file.
* `variance.txt`: Estimated variance of the LTQs *x<sub>gc</sub>* across cells *c* for each gene *g*. Note that these variances are different, and generally larger, than what one would obtain when directly calculating the variance of the estimates of *x<sub>gc</sub>* from the file `log_transcription_quotients.txt`. This is because the estimates in this file take into account the uncertainty on the estimates of the *x<sub>gc</sub>*. Thus, when estimates of true gene expression variability are needed, you are strongly advised to use the results in this file. Order corresponds to the order of genes in the `geneID.txt` file.
* `delta.txt`: Matrix of inferred log-fold changes *&delta;<sub>gc</sub> = x<sub>gc</sub>-&mu;<sub>g</sub>* for each gene *g* in each cell *c*. TSV file, columns correspond to cell IDs from `cellID.txt` file, rows correspond to gene IDs from `geneID.txt` file.
* `d_delta.txt`: Matrix of error-bars for the inferred log fold-changes *&delta;<sub>gc</sub>*. TSV file, columns correspond to cell IDs from `cellID.txt` file, rows correspond to gene IDs from `geneID.txt` file.
* `likelihood.txt`: This file encodes the posterior distribution of each gene's true variance in log-expression. For the numerical calculation of this distribution, the variance is a priori assumed to lie in the range *[v<sub>min</sub>,v<sub>max</sub>]*, and the distribution is evaluated on a grid of *N<sub>b</sub>* values spaced uniformly on a logarithmic scale. The file contains a matrix with likelihoods for each gene *g* and each grid point *b*. The first line lists the grid values *v<sub>b</sub> = v<sub>min</sub> exp(b &middot; &Delta;v)* for *b = 0,...,N<sub>b</sub>-1* with *&Delta;v = log(v<sub>max</sub>/v<sub>min</sub>)/(N<sub>b</sub>-1)*, so the first and last values are *v<sub>min</sub>* and *v<sub>max</sub>* themselves; the following lines correspond to the likelihood values of genes.

  | | | | | |
  |:-------|:------:|:------:|:------:|:---:|
  | Variance | 0.001000 | 0.001070 | 0.001146 | ... |
  | Gene 1 | 0.018 | 0.019 | 0.020 | ... |
  | Gene 2 | 0.0006 | 0.0051 | 0.0031 | ... |
  | ... | ... | ... | ... | ... |
  
## Usage
```
  ./Sanity <option(s)>
  Options:
    -h, --help                          Show this help message
    -v, --version                       Show the current version
    -f, --file                          Specify the input transcript count text file (.mtx for Matrix Market File Format)
    -mtx_genes, --mtx_gene_name_file    Specify the gene name text file (only needed if .mtx input file)
    -mtx_cells, --mtx_cell_name_file    Specify the cell name text file (only needed if .mtx input file)
    -d, --destination                   Specify the destination path (default: pwd)
    -n, --n_threads                     Specify the number of threads to be used (default: 4)
    -e, --extended_output               Option to print extended output (default: false, choice: false,0,true,1)
    -v_m, --v_method                    Option to specify the method for variance estimation (default: MAP, choice: MAP, EAP, MLE, MARG)
    -vmin, --variance_min               Minimal value of variance in log transcription quotient (default: 0.001)
    -vmax, --variance_max               Maximal value of variance in log transcription quotient (default: 50)
    -nbin, --number_of_variance_bins    Number of bins for the variance in log transcription quotient  (default: 160)
    -no_norm, --no_cell_size_normalization  Option to skip cell size normalization (default: false, choice: false,0,true,1)
```

## Installation
* Clone the GitHub repository
```
git clone https://github.com/jmbreda/Sanity.git
```
* Install OpenMP library
	* On Linux  
	If not already installed (Check with `ldconfig -p | grep libgomp`, no output if not installed), do
	```
	sudo apt-get update
	sudo apt-get install libgomp1
	```
	* On mac OS  
	Apple's own `g++`/`clang++` does not support OpenMP, so you need a real GCC.

		* Using brew (recommended)  
		```
		brew install gcc
		```
		The Makefile detects this automatically and picks the right `g++` binary, no further changes needed.

		* Using macports  
		```
		port install gcc13
		```
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tell the Makefile to use it by compiling with `make CC=g++-mp-13` (adjust the version number to the package you installed).
* Move to the source code directory and compile.
```
cd Sanity/src
make
```
* The binary file is located in
```
Sanity/bin/Sanity
```
* Alternatively, pre-built binaries are published on the repository's [Releases](https://github.com/jmbreda/Sanity/releases) page.

## Sanity_distance
Compute cell-cell distances from Sanity output files. Needs extended outputs of Sanity (`-e 1` option).
### Input
* `-f`: The output folder of the Sanity run, specified with the `-d` option in Sanity (`'path/to/folder'`)
* `-s2n`: (optional) The gene signal to noise ratio used as gene cut-off (double, default: `1.0`)
* `-err`: (optional) Compute cell-cell distance with or without error bars (boolean, default: `1` or `true`)
* `-n`: (optional) Number of threads (integer, default: `4`)
### Output
* Cell-cell distance: *(N<sub>c</sub>(N<sub>c</sub>-1)/2)* vector of cell to cell distances *dist(cell<sub>i</sub>,cell<sub>j</sub>), i=1,...,N<sub>c</sub>-1, j=i+1,...,N<sub>c</sub>*, with *N<sub>c</sub>* the number of cells. 
  ||
  |:------:|
  |*dist(cell<sub>1</sub>,cell<sub>2</sub>)*|
  |*dist(cell<sub>1</sub>,cell<sub>3</sub>)*|
  |*dist(cell<sub>1</sub>,cell<sub>4</sub>)*|
  |...|
  |*dist(cell<sub>N<sub>c</sub>-2</sub>,cell<sub>N<sub>c</sub>-1</sub>)*|
  |*dist(cell<sub>N<sub>c</sub>-2</sub>,cell<sub>N<sub>c</sub></sub>)*|
  |*dist(cell<sub>N<sub>c</sub>-1</sub>,cell<sub>N<sub>c</sub></sub>)*|
located in the Sanity output folder (specified with `-f` option), named `cell_cell_distance_[...].txt`, depending on the `-err` and `-s2n` options.
### Usage
```
./Sanity_distance <option(s)>
Options:
	-h,--help		Show this help message
	-v,--version		Show the current version
	-f,--folder		Specify the input folder with extended output from Sanity
	-s2n,--signal_to_noise_cutoff	Minimal signal/noise of genes to include in the distance calculation (default: 1.0)
	-err,--with_error_bars	Compute cell-cell distance taking the errorbar epsilon into account (default: true)
	-n,--n_threads		Specify the number of threads to be used (default: 4)
```

### Installation
Same dependencies as Sanity (see above).

* Move to the source code directory and compile.
```
cd Sanity/src
make Sanity_distance
```
* The binary file is located in
```
Sanity/bin/Sanity_distance
```

<a name="sorting-mtx-file-by-row-gene-indices"></a>
### Sorting MTX file by row (gene) indices

Use script `sort_mtx_by_row.py` in the `Sanity/scripts` folder to sort the MTX file by row (gene) indices. The script takes as input the MTX file. The output is a sorted MTX file by default saved as `sorted_[original filename].gz` in the same folder as the input file. User can specify the output file name and path with the `-o` option. The saved MTX file is compressed with gzip.

Requirements: Python3, sort

Example usage:

```bash
python3 sort_mtx_by_row.py -i input.mtx
```

## Tests
A simple test comparing Sanity's output across the `-v_m` methods (MAP, EAP, MLE, MARG) on a small count matrix is available in the `tests` directory. See [tests/README.md](tests/README.md) for details on running it.

## Help
For any questions or assistance regarding Sanity, please post your question in the issues section.
