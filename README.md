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
      

## Calculating ERCs on the command line
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


## Calculating ERCs in Python
The pipeline is easily accessible in Python. For practical examples of usage, check out the source of `cli.py`.

Setting up the environment:
```python
from pipeline import ErcWorkspace

# The working dir and tree args are the only required arguments
workspace = ErcWorkspace("directory/to/run/in", "path/to/tree/topology.nwk")  
```
Note that this immediately changes the Python's runtime current working directory to the working directory passed as the
first argument.

Definitions of the additional optional arguments;
* segmented: If True, run ERC on pieces of the alignments
* segment_size: The size of the pieces
* internal_requirement: The minimum required size of internal branches to perform correlations with (only applies to 
  normal ERCs and by default only all terminal branches are included in ERC calculations)
* include_terminal: If you pass an internal_requirement argument, should terminal branches also be included in ERC 
  calculations?
* recalculate: If True, even if ERCs have been previously calculated for tree pairs in the environment, recalculate the 
  correlation results.
* sliding_window: If True, ERCs run on pieces of alignments are normalized with a sliding window.
* skip_align: If True, skip alignment of the input sequences.
* skip_trim: If True, skip trimming of alignments.
* taxon_set: Pass a list of taxa to perform calculations with, by default all taxa are considered.
* time_corrected: If True, instead of standard Spearman's test, run Spearman's partial correlations controlling for time.
* id2name: A dictionary that if passed, will be used to convert protein ids to readable names in outputs.

Include previously run ERCs (assuming you ran the pipeline in directory: `old/erc_dir`):
```python
from pipeline import register_erc_datasource

register_erc_datasource("old/erc_dir/tree/", "old/erc_dir/ercs.csv")
```

Add alignments to the workspace:
```python
# Add alignments for all-by-all calculations
workspace.add_alignment("path/to/alignment.fasta")

# Add specific pairs of alignments to run calculations on
workspace.add_alignment("path/to/alignment1.fasta", "path/to/alignment2.fasta")

# Add alignments that should be concatenated together when run for ERCs
workspace.add_concatenated_alignment(["path/to/alignment1.fasta", "path/to/alignment2.fasta"], "concatenated_name")

# You can specify specific pairs of concatenated pairs as well
workspace.add_concatenated_alignment(["path/to/alignment1.fasta", "path/to/alignment2.fasta"], "concatenated_name",
                                     ["path/to/alignment3.fasta", "path/to/alignment4.fasta"], "concatenated_name2")
```

Run the ERCs:
```python
# If running outside of a coroutine:
import asyncio
asyncio.get_event_loop().run_until_complete(workspace.run())

# If inside a coroutine:
await workspace.run()
```

Following calculations, you can load the data into a `networkx` `Graph` object using:
```python
net = workspace.generate_network()
```
This network has each protein ran as nodes, with edges connecting the nodes based on ERCs. Each
edge has the `rho` and `p` properties representing the ERC results.

Save data:
```python
# Save the data to an excel spreadsheet
workspace.export_results("filename.xlsx")
```

Access Mammalian TimeTree:
```python
from pipeline import _10mya_cutoff, _20mya_cutoff, _30mya_cutoff, prune_tree
from utilities import _self_path, safe_phylo_read
import os.path as osp

tree = safe_phylo_read("path/to/file.nwk")  # Read general newick trees
timetree = safe_phylo_read(osp.join(_self_path(), 'data', 'finished_mam_timetree.nwk'))  # Read the mammalian TimeTree

# Prune mammalian tree to the 10my, 20my, and 30my taxa sets respectively
_10my = prune_tree(timetree, _10mya_cutoff)
_20my = prune_tree(timetree, _20mya_cutoff)
_30my = prune_tree(timetree, _30mya_cutoff)
```

### Misc. Features in Python
Improved async performance:
```python
# If the uvloop package is installed (unix only), you can run the following to improve performance before running ERCs:
from utilities import try_hook_uvloop

try_hook_uvloop()
```

Extract rate data from trees:
```python
from pipeline import get_rates

# Calculate rates for the trees passed. First argument is the PhyloTree object for the species topology.
# Second arg determines whether to prune the tree topologies so all trees have the same topologies
# Third arg is a list of taxa to limit calculations to
# The remaining args are the trees for proteins of interest
taxa, rates = await get_rates(tree_topology_object, True, None, tree1, tree2)
# Taxa is a list of taxa names
# Rates is a list of lists. The first list are all the time units, the following lists are the rates. 
# With the indices of each element corresponding to the taxon in the matching index of the taxa list.
```

Convert rate data to correlations:
```python
from pipeline import get_rates, rates_to_correlation

rate_info = await get_rates(tree_topology_object, True, None, tree1, tree2)
rho, p = await rates_to_correlation(rate_info)
```

Run an enrichment analysis on protein sets:
```python
from pipeline import enrich_network

# Where protein_symbols is a list of protein identifiers (ids or protein symbols),
# background_protein_symbols is a list of protein identifiers from the background set (ids or protein symbols),
# id2name is a dictionary (can be empty) mapping ids to symbols, and the last argument is the base file name for 
# enrichment reports.
enrich_network(protein_symbols, background_protein_symbols, id2name, "enrichment_results")
```

*The following features assume you generated a network object from your ERCs*

Generate Reciprocal-Rank 20 Networks (RRNs)
```python
from utilities import rrn_nets

# Each step represents an intermediate RRN result
step1, step2, step3 = rrn_nets(net, "central_protein")  # You can pass multiple central proteins
```

Export a network diagram (requires pygraphviz):
```python
from utilities import graphviz_network_plot

# You can pass the following optional arguments:
# * highlight: A dictionary mapping proteins -> color code for nodes
# * circo: If True, use the circo graphviz layout instead of the default
# * highlight_by_font: If True, highlighted nodes will have the text highlighted instead of the background
graphviz_network_plot(my_network, "my_network.png")
```

## Citations

Please cite the ERC method with the following reference:
```
Yan, Zhichao, Gongyin Ye, and John H. Werren. "Evolutionary rate correlation between mitochondrial-encoded and mitochondria-associated nuclear-encoded proteins in insects." Molecular biology and evolution 36.5 (2019): 1022-1036.
```

Please cite this specific implementation of ERCs with the following reference:
```
TODO
```
