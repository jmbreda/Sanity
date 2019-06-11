# Sanity

**Sa**mpling **N**oise based **I**nference of **T**ranscription Activit**Y** : Filtering of Poison noise on a single-cell RNA-seq UMI count matrix

Sanity infers the log expression levels *x<sub>gc</sub>* of gene *g* in cell *c* by filtering out 
the Poisson noise on the UMI count matrix *n<sub>gc</sub>* of gene *g* in cell *c*.

## Input

* Count Matrix : *(N<sub>g</sub> x N<sub>c</sub>)* matrix with *N<sub>g</sub>* the number of genes and *N<sub>c</sub>* the number of cells. tsv or csv format.

| GeneID | Cell 1 | Cell 2 | Cell 3 | ...
|:-------|:------:|:------:|:------:|------:|
| Gene 1 | 1.0 | 3.0 | 0.0 |
| Gene 2 | 2.0 | 6.0 | 1.0 |
| ... | |

* (optional) Destination folder (`'path/to/output/folder'`).
* (optional) Number of threads (integer number).
* (optional) Print extended output (Boolean, `'true', 'false', '1'` or `'0'`)

## Output

* expression_level.txt : *(N<sub>g</sub> x N<sub>c</sub>)* table of inferred log expression levels

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.25 | -0.29 | -0.54 |
  | Gene 2 | -0.045 | -0.065 | 0.11 |
  | ... | |

* d_expression_level.txt : *(N<sub>g</sub> x N<sub>c</sub>)* table of error bars on inferred log expression levels

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
	-e,--extended_output	Option t print extended output (default: false)
```

## Installation
* clone the GitHub repository
```
git clone https://github.com/jmbreda/Sanity.git
```
* intall needed library
```
sudo apt-get update
sudo apt-get install libgomp1
```
* Move to the source code directory
```
cd Sanity/src
```
* Compile the code.
```
make
```
* The binary file is located in
```
Sanity/bin/Sanity
```
