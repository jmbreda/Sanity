# Sanity

**Sa**mpling **N**oise based **I**nference of **T**ranscription Activit**Y** : Filtering of Poison noise on a single-cell RNA-seq UMI count matrix

Single-cell RNA sequencing normalization algorithm presented in the publication [Bayesian inference of gene expression states from single-cell RNA-seq data -
J Breda, M Zavolan, E van Nimwegen - Nature Biotechnology, 2021](https://www.nature.com/articles/s41587-021-00875-x).

Sanity infers the log expression levels *x<sub>gc</sub>* of gene *g* in cell *c* by filtering out 
the Poisson noise on the UMI count matrix *n<sub>gc</sub>* of gene *g* in cell *c*.



### Reproducibility
The raw UMI count and normalized datasets mentioned in benchmarking in the associated [publication](https://www.nature.com/articles/s41587-021-00875-x) are available on [![DO I](https://zenodo.org/badge/DOI/10.5281/zenodo.4009187.svg)](https://zenodo.org/record/4009187). Files are named [*dataset name*]\_UMI\_counts.txt.gz and [*dataset name*]\_[*tool name*]\_normalization.txt.gz.

The scripts used for running the bechmarked normalization methods and for making the figures of the preprint are in the reproducibility folder.

## Input

* UMI count matrix: *(N<sub>g</sub> x N<sub>c</sub>)* matrix with *N<sub>g</sub>* the number of genes and *N<sub>c</sub>* the number of cells. Format: tab-separated, comma-separated, or space-separated values. (`'path/to/text_file'`)

| GeneID | Cell 1 | Cell 2 | Cell 3 | ...
|:-------|:------:|:------:|:------:|------:|
| Gene 1 | 1.0 | 2.0 | 0.0 |
| Gene 2 | 6.0 | 3.0 | 1.0 |
| ... | |

* (Alternatively) Matrix Market File Format: Sparse matrix of UMI counts. Automatically recognized by `.mtx` extension of the input file. Named `matrix.mtx` by cellranger 2.1.0 and 3.1.0 (10x Genomics). (`'path/to/text_file.mtx'`)
	* (optional) Gene ID file: Named `genes.tsv` by cellranger 2.1.0 and `features.tsv` by cellranger 3.1.0 (10x Genomics). (`'path/to/text_file'`)
	* (optional) Cell ID file: Named `barcodes.tsv` by cellranger 2.1.0 and 3.1.0 (10x Genomics).  (`'path/to/text_file'`)
* (optional) Destination folder (`'path/to/output/folder'`, default: `pwd`)
* (optional) Number of threads (integer, default: `4`)
* (optional) Print extended output (Boolean, `'true', 'false', '1'` or `'0'`, default: `false`)
* (optional) Print additional output (with suffix "_vmax.txt") obtained using the maximum posterior estimate for the gene-variances $v_g$, rather than integrating over the posterior of $v_g$. This is the most correct output when one wants to reconstruct the likelihood expression from the posterior estimates, for example in the distance-calculation script. (string, `'true', 'false', '1', '0'` or `'only_max_output'`, default: `false`)
* (optional) Minimal and maximal considered values of the variance in log transcription quotients (double, default: *v<sub>min</sub>=*`0.001` *v<sub>max</sub>=*`50`)
* (optional) Number of bins for the variance in log transcription quotients (integer, default: `160`)
* (optional) Option to skip cell size normalization (Boolean, `'true', 'false', '1'` or `'0'`, default: `false`)
## Output

* log_transcription_quotients.txt: This file contains the estimated values of the log-transcription quotients (LTQs) for each gene in each cell. The LTQ *x<sub>gc</sub>* of gene *g* in cell *c* corresponds to the estimated logarithm of the fraction of mRNAs in cell *c* that belong to gene *g*. The LTQs are thus normalized such that *&Sigma;<sub>g</sub> exp(x<sub>gc</sub>) = 1* for each cell *c*. In order to get an estimate of the number of mRNAs for gene *g* in cell *c* one would thus need to multiply *exp(x<sub>gc</sub>)* by the estimated total number of mRNAs *M* in the cell.

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | -13.7227 | -13.722 | -13.729 |
  | Gene 2 |  -9.96744 | -10.2522 | -10.1453 |
  | ... | |
  
* ltq_error_bars.txt : Table with the error-bars on the estimates of the LTQs *x<sub>gc</sub>* for each gene *g* in each cell *c*.

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.630111 | 0.630198 | 0.624802 |
  | Gene 2 | 0.315551 | 0.325912 | 0.301861 |
  | ... | |

* (optional) log_transcription_quotients_vmax.txt: Analogous file to log_transcription_quotients.txt but reporting the LTQs obtained by using the maximum-posterior gene-variance $v_g$. This file is only produced when `-max v` is set to `true` or `only_max_output`.
* (optional) ltq_error_bars_vmax.txt: See the description for log_transcription_quotients_vmax.txt.

## Extended output (optional)

* mu.txt : Estimated average LTQ *&mu;<sub>g</sub>* of each gene *g* (averaged over all cells)
* d_mu.txt : Error bars on the inferred mean LTQs *&mu;<sub>g</sub>*.
* variance.txt : Estimated variance of the LTQs *x<sub>gc</sub>* across cells *c* for each gene *g*. Note that these variances are different, and generally larger, than what one would obtain when directly calculating the variance of the estimates of *x<sub>gc</sub>* from the file log_transcription_quotients.txt. This is because the estimates in this file take into account the uncertainty on the estimates of the *x<sub>gc</sub>*. Thus, when estimates of true gene expression variability are needed, you are strongly adviced to use the results in this file.
* delta.txt : Matrix of inferred log-fold changes *&delta;<sub>gc</sub> = x<sub>gc</sub>-&mu;<sub>g</sub>* for each gene *g* in each cell *c*.
* d_delta.txt : Matrix of error-bars for the inferred log fold-changes *&delta;<sub>gc</sub>*.
* ..._vmax.txt: All the above files will be obtained with a `_vmax.txt`-suffix when `-max v` is set to `true` or `only_max_output`.
* likelihood.txt : This file encodes the posterior distribution of each gene’s true variance in log-expression. For the numerical calculation of this distribution, the variance is a prior assumed to lie in the range *[v<sub>min</sub>,v<sub>max</sub>]* and is discretized into *N<sub>b</sub>* bins uniformly on a logarithmic scale. The file contains the matrix with posterior values *P<sub>gb</sub>* for each gene *g* and each bin *b*.

  | | | | | |
  |:-------|:------:|:------:|:------:|------:|
  |Variance | 0.01 | 0.0107 | 0.0114 | ... |
  | Gene 1 | 0.018 | 0.019 | 0.020 |
  | Gene 2 | 0.0006 | 0.0051 | 0.0031 |
  |...|
  
## Usage
```
  ./Sanity <option(s)> SOURCES
  Options:
	-h,--help		Show this help message
	-v,--version		Show the current version
	-f,--file		Specify the input transcript count text file (.mtx for Matrix Market File Format)
	-mtx_genes,--mtx_gene_name_file	Specify the gene name text file (only needed if .mtx input file)
        -mtx_cells,--mtx_cell_name_file	Specify the cell name text file (only needed if .mtx input file)
	-d,--destination	Specify the destination path (default: pwd)
	-n,--n_threads		Specify the number of threads to be used (default: 4)
	-e,--extended_output	Option to print extended output (default: false, choice: false,0,true,1)
	-max_v,--get_output_for_maxlik_variance	Option to obtain additional output obtained for maximum-posterior gene-variance (default: false, choice: false,0,true,1,only_max_output)
	-vmin,--variance_min	Minimal value of variance in log transcription quotient (default: 0.001)
	-vmax,--variance_max	Maximal value of variance in log transcription quotient (default: 50)
	-nbin,--number_of_bins	Number of bins for the variance in log transcription quotient  (default: 160)
	-no_norm,--no_cell_size_normalization	Option to skip cell size normalization (default: false, choice: false,0,true,1)
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
	
	* On mac OS using macports  
	Install the `gcc9` package
	```
	port install gcc9
	```
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Change the first line of `src/Makefile` from `CC=g++` to `CC=g++-mp-9`
	
	* On mac OS using brew  
	Install the `gcc9` package  
	```
	brew install gcc9
	```
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Change the first line of `src/Makefile` from `CC=g++` to `CC=g++-9`
* Move to the source code directory and compile.
```
cd Sanity/src
make
```
* The binary file is located in
```
Sanity/bin/Sanity
```
* Alternatively, the already compiled binary for macOS is located in
```
Sanity/bin/Sanity_macOS
```

## Sanity_distance
Compute cell-cell distances from Sanity output files. Needs extended outputs of Sanity (`-e 1` option).
### Input
* The output folder of the Sanity run, specifiied with the `-d` option in Sanity (`'path/to/folder'`)
* (optional) The gene signal to noise ratio used as gene cut-off (double, default: `1.0`)
* (optional) Compute distances with or without errorbars (boolean, default: `1` or `true`)
* (optional) Number of threads (integer, default: `4`)
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
./Sanity_distance <option(s)> SOURCES
Options:
	-h,--help		Show this help message
	-v,--version		Show the current version
	-f,--folder		Specify the input folder with extended output from Sanity
	-s2n,--signal_to_noise_cutoff	Minimal signal/noise of genes to include in the distance calculation (default: 1.0)
	-err,--with_error_bars	Compute cell-cell distance taking the errobar epsilon into account (default: true)
	-n,--n_threads		Specify the number of threads to be used (default: 4)
```
### Instalation
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
## Help
For any questions or assistance regarding Sanity, please post your question in the issues section.
