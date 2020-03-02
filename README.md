[![Build Status](https://travis-ci.org/tseemann/cgmlst-dists.svg?branch=master)](https://travis-ci.org/tseemann/cgmlst-dists)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: C99](https://img.shields.io/badge/Language-C99-orangered.svg)](https://en.wikipedia.org/wiki/C99)
<!-- ![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.1411986.svg) -->

# cgmlst-dists

Calculate distance matrix from cgMLST allele call tables of ChewBBACA

## Quick Start

```
% cat test/boring.tab

FILE    G1      G2      G3      G4      G5      G6
S1      1       2       3       2       1       5
S2      1       1       1       1       NIPH    5
S3      1       2       3       4       1       3
S4      1       LNF     2       4       1       3
S5      1       2       ASM     2       1       3

% cgmlst-dists test/boring.tab > distances.tab

This is cgmlst-dists 0.1.0
Loaded 5 samples x 6 allele calls
Calulating distances... 100.00%
Done.

% cat distances.tab

        S1      S2      S3      S4      S5
S1      0       3       2       3       1
S2      3       0       4       3       3
S3      2       4       0       1       1
S4      3       3       1       0       1
S5      1       3       1       1       0
```

## Installation

`cgmlst-dists` is written in C and has no other dependencies.

### Homebrew
```
brew install brewsci/bio/cgmlst-dists  # COMING IN APRIL 2020
```

### Bioconda
```
conda install -c bioconda -c conda-forge cgmlst-dists  # COMING SOON
```

### Source

```
git clone https://github.com/tseemann/cgmlst-dists.git
cd cgmlst-dists
make

# run tests
make check

# optionally install to a specific location (default: /usr/local)
make PREFIX=/usr/local install
```

## Options

### `cgmlst-dists -h` (help)

```
SYNOPSIS
  Pairwise CG-MLST distance matrix from allele call tables
USAGE
  cgmlst-dists [options] chewbbaca.tab > distances.tsv
OPTIONS
  -h    Show this help
  -v    Print version and exit
  -q    Quiet mode; do not print progress information
  -c    Use comma instead of tab in output
URL
  https://github.com/tseemann/cgmlst-dists
```

### `cgmlst-dists -v` (version)

Prints the name and version separated by a space in standard Unix fashion.

```
cgmlst-dists 0.1.0
```

### `cgmlst-dists -q` (quiet mode)

Don't print informational messages, only errors.

### `cgmlst-dists -c` (CSV mode)

Use a comma instead of a tab in the output table.



## Issues

Report bugs and give suggesions on the
[Issues page](https://github.com/tseemann/cgmlst-dists/issues)

## Related software

* [chewBBACA](https://github.com/B-UMMI/chewBBACA)
* [snp-dists](https://github.com/tseemann/snp-dists)

## Licence

[GPL Version 3](https://raw.githubusercontent.com/tseemann/cgmlst-dists/master/LICENSE)

## Authors

* [Torsten Seemann](https://github.com/tseemann)
