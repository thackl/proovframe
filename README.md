[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/proovframe/README.html)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/proovframe.svg?style=flat)](https://anaconda.org/bioconda/proovframe)
[![DOI](https://img.shields.io/badge/DOI-10.1101/2021.08.23.457338-blue)](https://doi.org/10.1101/2021.08.23.457338)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/proovframe/badges/license.svg)](https://anaconda.org/bioconda/proovframe)


proovframe: frame-shift correction for long read (meta)genomics
=========================================

Gene prediction on long reads, aka PacBio and Nanopore, is often impaired by
indels causing frameshift. Proovframe detects and corrects frameshifts in coding
sequences from raw long reads or long-read derived assemblies.  

<img src="implementation.png" width="1000px" height="600px" />

Proovframe uses frameshift-aware alignments to reference proteins as guides, and
conservatively restores frame-fidelity by 1/2-base deletions or insertions of
`N/NN`s, and masking of premature stops (`NNN`).

Good results can already be obtained with distantly related guide proteins-
successfully tested with sets with <60% amino acid identity.

It can be used as an additional polishing step on top of classic
consensus-polishing approaches for assemblies.

It can be used on raw reads directly, which means it can be used on data lacking
sequencing depth for consensus polishing - a common problem for a lot of rare
things from environmental metagenomic samples, for example.

## Install 

### bioconda

```
conda install -c bioconda proovframe
```

### Manual

Requires  [DIAMOND v2.0.3](https://github.com/bbuchfink/diamond) or newer for mapping.

```
git clone https://github.com/thackl/proovframe
```

It is ready to be used. The tool lives in `proovframe/bin/proovframe`


## Usage

map proteins to reads: 
``` 
proovframe/bin/proovframe map -a proteins.faa -o raw-seqs.tsv raw-seqs.fa
```

fix frameshifts in reads:
```
proovframe/bin/proovframe fix -o corrected-seqs.fa raw-seqs.fa raw-seqs.tsv
```


## Citing

If you use proovframe and DIAMOND please cite: 

- Hackl T, Trigodet F, Murat Eren A, Biller SJ, Eppley JM, Luo E, et al. "proovframe: frameshift-correction for long-read (meta)genomics", bioRxiv. 2021. p. 2021.08.23.457338. doi:10.1101/2021.08.23.457338
- Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
