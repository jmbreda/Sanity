# TrueVar

TrueVar : Filtering of Poison noise on a single-cell RNA-seq UMI count matrix

TrueVar infers the log expression levels x_gc of gene g in cell c by filtering out 
the Poisson noise on the UMI count matrix n_gc of gene g in cell c.

## Input

  *Count Matrix* : *(N_g x N_c)* matrix with *N_g* the number of genes and *N_c* the number of cells

| GeneID | Cell 1 | Cell 2 | Cell 3 | ...
|:-------|:------:|:------:|:------:|------:|
| Gene 1 | 1.0 | 3.0 | 0.0 |
| Gene 2 | 2.0 | 6.0 | 1.0 |
| ... | |

## Output

* *expression_level.txt* : *(N_g x N_c)* table of infered log expression levels

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.25 | -0.29 | -0.54 |
  | Gene 2 | -0.045 | -0.065 | 0.11 |
  | ... | |

* d_expression_level.txt : *(N_g x N_c)* table of error bars on infered log expression levels

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ...
  |:-------|:------:|:------:|:------:|------:|
  | Gene 1 | 0.015 | 0.029 | 0.042 |
  | Gene 2 | 0.0004 | 0.0051 | 0.0031 |
  | ... | |

## Extended output (optional)

* mu.txt : *(N_g x 1)* vector of infered mean log expression levels
* d_mu.txt : *(N_g x 1)* vector of infered error bars on mean log expression levels
* delta.txt : *(N_g x N_c)* matrix of infered log expression levels centered in 0
* d_delta.txt : *(N_g x N_c)* matrix of infered error bars log expression levels centered in 0
* likelihood.txt : *(N_g+1 x N_b)* matrix of normalized variance likelihood per gene, with *N_b* the number of bins on the variance.

  | | | | | |
  |:-------|:------:|:------:|:------:|------:|
  |Variance | 0.01 | 0.0107 | 0.0114 | ... |
  | Gene 1 | 0.018 | 0.019 | 0.020 |
  | Gene 2 | 0.0006 | 0.0051 | 0.0031 |
  |...|
  
## Usage
```
  ./TrueVar <option(s)> SOURCES
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
git clone https://github.com/jmbreda/TrueVar.git
```
* intall needed library
```
sudo apt-get update
sudo apt-get install libgomp1
```
* Move to the source code directory
```
cd TrueVar/src
```
* Compile the code.
```make```
* The binary is file is located in ```TrueVar/bin```
