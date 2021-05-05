# ERC-Pipeline

The pipeline in this repository allows for the generation of evolutionary rate correlations (ERCs). ERCs use a 
phylogeny-based approach to look for evolutionary signatures of potential protein-protein interaction. Since ERCs can
be fairly straightforward to produce, they can be used as an easy method of discovering candidate protein-protein 
interactions.

## Installation
### Environment Requirements
The pipeline should be compatible with Python version 3.7 or higher. 
R scripts were written for R 3.6.1 ("Action of the Toes") but are optional.
Note that the pipeline was made to take advantage of 
multiple cores. We recommend at least 5 available cores and 16GB of RAM for large data sets.

### Python Package Dependencies
The external Python dependencies can be installed using pip (most packages can also be found on conda/bioconda):
`pip install -r requirements.txt`

Note that `uvloop` and `pygraphviz` are *optional* dependencies so if they fail to install on your system, you can still
run the pipeline. Just note that `pygraphviz` is required for network diagrams.

### R Package Requirements
All R packages are optional. They can be installed using the `install.packages()` function in base R.

* If using time-corrected partial correlation calculations, `ppcor`
* If using any of the `generate_figs.R` script code, `ggplot2, readxl, writexl, ggvenn, dplyr, gridExtra, ape, picante, EnvStats, Cairo (on Windows)` 

### External Tooling Requirements
All of the following are *required* to be installed.
* MAFFT (https://mafft.cbrc.jp/alignment/software/) - For generating alignments.
* trimAl (http://trimal.cgenomics.org/) - For trimming alignments.
* IQ-TREE 2 (http://www.iqtree.org/) - For calculating phylogenetic trees.

### Preparation of Data
1. Time-scaled phylogeny (Newick format)
    * If not specified, the pipeline will use the Mammalian time-scaled phylogeny generated using TimeTree 
      (http://timetree.org/) for all the Mammalian taxa available in OrthoDB v10 (https://www.orthodb.org/). 
      **Implementation Note:** The species names in this case are formatted in all caps with spaces replaced with 
      underscores (i.e. `Homo sapiens` is formatted as `HOMO_SAPIENS`). You can find the full file in 
      [data/finished_mam_timetree.nwk](./data/finished_mam_timetree.nwk).
2. Protein sequences (FASTA format)
    * All proteins of interest should be provided as FASTA formatted sequences (these are assumed to be unaligned) in a
      directory. Note that the title of each sequence must be the taxon name for each sequence corresponding to the 
      time-scaled phylogeny (i.e. if you have a human protein sequence, the FASTA sequence title must be exactly 
      `HOMO_SAPIENS` if using the default tree).
    * This pipeline does not attempt to disambiguate paralogous sequences. So each protein sequence should be singe-copy,
      any multi-copy sequences in the FASTA files will be totally ignored.
      

## Calculating ERCs
Assuming the sequences are in a directory called `ALIGNMENT_DIRECTORY` and the time-scaled phylogeny is a file called
`PHYLOGENY.nwk` and you wish to run calculations in a directory called `OUTPUT_DIRECTORY`.

`python3 cli.py --timetree PHYLOGENY.nwk --alignments ALIGNMENT_DIRECTORY --wd OUTPUT_DIRECTORY`

* Note that if you are using the default mammalian phylogeny, you can omit the `--timetree PHYLOGENY.nwk` argument. 
  Additionally, you would be able to use the 20MY or 30MY ERC calculations described in Varela et al 2021 using either
  `--erc-type 20my` or `--erc-type 30my` respectively. 
  
* If you wish to use the time-corrected correlation method, you can add the `--erc-type bt` argument (requires R).

* If you want to only run ERCs on specific pairs of proteins, you can pass the 
  `--align-pair /path/to/align1.fasta /path/to/align2.fasta` option as many times as needed.

* If your sequence data are already prepared, you can use `--skip-align` and `--skip-trim` to skip alignment and 
  trimming of alignments, respectively.
  
* If your alignments are titled based on non-readable identifiers instead of protein names/symbols (for example: 
  OrthoDB ids), you can create a tab-separated file with 2 columns: "alignment_identifier" and "readable_name". This
  will be used to replace the original alignment names with the readable names if you pass the path to this file to the
  argument: `--id2name /path/to/file.tsv`

* If you want to run ERCs along pieces of the alignments, add the `--segment` argument. This can be modified using the 
  `--slide` argument, which rather than splitting the alignment into kmers, will normalize the data using a sliding 
  window based on kmers. You can also use `--kmer K` to change the size of the kmers (replace "K" with the number). 


## Citations

Please cite the ERC method with the following reference:
```
Yan, Zhichao, Gongyin Ye, and John H. Werren. "Evolutionary rate correlation between mitochondrial-encoded and mitochondria-associated nuclear-encoded proteins in insects." Molecular biology and evolution 36.5 (2019): 1022-1036.
```

Please cite this specific implementation of ERCs with the following reference:
```
TODO
```
