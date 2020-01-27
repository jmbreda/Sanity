# Sanity

**Sa**mpling **N**oise based **I**nference of **T**ranscription Activit**Y** : Filtering of Poison noise on a single-cell RNA-seq UMI count matrix

Sanity infers the log expression levels *x<sub>gc</sub>* of gene *g* in cell *c* by filtering out 
the Poisson noise on the UMI count matrix *n<sub>gc</sub>* of gene *g* in cell *c*.

See our [preprint](https://doi.org/10.1101/2019.12.28.889956 "bioRxiv: Bayesian inference of the gene expression states of single cells from scRNA-seq data") for more details.

### Reproducibility
The raw and normalized datasets mentionned in the [preprint](https://doi.org/10.1101/2019.12.28.889956 "bioRxiv: Bayesian inference of the gene expression states of single cells from scRNA-seq data") are available on [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3597622.svg)](https://doi.org/10.5281/zenodo.3597622). Files are named [*dataset name*]\_UMI\_counts.txt.gz and [*dataset name*]\_[*tool name*]\_normalization.txt.gz.

The scripts used for running the bechmarked normalization methods and for making the figures of the preprint are in the reproducibility folder.

## Input

* Count Matrix : *(N<sub>g</sub> x N<sub>c</sub>)* matrix with *N<sub>g</sub>* the number of genes and *N<sub>c</sub>* the number of cells. tsv or csv format.

| GeneID | Cell 1 | Cell 2 | Cell 3 | ...
|:-------|:------:|:------:|:------:|------:|
| Gene 1 | 1.0 | 3.0 | 0.0 |
| Gene 2 | 2.0 | 6.0 | 1.0 |
| ... | |

* (optional) Destination folder (`'path/to/output/folder'`)
* (optional) Number of threads (integer)
* (optional) Print extended output (Boolean, `'true', 'false', '1'` or `'0'`)
* (optional) Minimal and maximal considered values of the variance in log transcription quotients (double)
* (optional) Number of bins for the variance in log transcription quotients (integer)
## Output

* log_transcription_quotients.txt : *(N<sub>g</sub> x N<sub>c</sub>)* table of inferred log expression levels. The gene expression levels are normalized to 1, meaning that the summed expression of all genes in a cell is approximately 1. Note that we use the natural logarithm, so to change the normalization one should multiply the exponential of the expression by the wanted normalization (*e.g.* mean or median number of captured gene per cell).

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.25 | -0.29 | -0.54 |
  | Gene 2 | -0.045 | -0.065 | 0.11 |
  | ... | |

* ltq_error_bars.txt : *(N<sub>g</sub> x N<sub>c</sub>)* table of error bars on inferred log expression levels

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.015 | 0.029 | 0.042 |
  | Gene 2 | 0.0004 | 0.0051 | 0.0031 |
  | ... | |

## Extended output (optional)

* mu.txt : *(N<sub>g</sub> x 1)* vector of inferred mean log expression levels
* d_mu.txt : *(N<sub>g</sub> x 1)* vector of inferred error bars on mean log expression levels
* variance.txt : *(N<sub>g</sub> x 1)* vector of inferred variance per gene in log expression levels
* delta.txt : *(N<sub>g</sub> x N<sub>c</sub>)* matrix of inferred log expression levels centered in 0
* d_delta.txt : *(N<sub>g</sub> x N<sub>c</sub>)* matrix of inferred error bars log expression levels centered in 0
* likelihood.txt : *(N<sub>g</sub>+1 x N<sub>b</sub>)* matrix of normalized variance likelihood per gene, with *N<sub>b</sub>* the number of bins on the variance.

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
	-f,--file		Specify the input transcript count text file
	-d,--destination	Specify the destination path (default: pwd)
	-n,--n_threads		Specify the number of threads to be used (default: 4)
	-e,--extended_output	Option to print extended output (default: false)
	-vmin,--variance_min	Minimal value of variance in log transcription quotient (default: 0.01)
	-vmax,--variance_max	Maximal value of variance in log transcription quotient (default: 20)
	-nbin,--number_of_bins	Number of bins for the variance in log transcription quotient  (default: 116)
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
* Move to the source code directory and compile. The compilation takes a few seconds.
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

## Help
For any questions or assistance regarding Sanity, please post your question the issues section or contact us at jeremie.breda@unibas.ch
