import asyncio
import itertools
import math
import os
import shutil
import sys
import traceback
from collections import namedtuple
import os.path as osp
from glob import glob
from typing import Union, List, Optional, Tuple, Dict

from ete3 import PhyloTree
import networkx as nx
import numpy as np
from numpy import array, nan
import gseapy as gp
from scipy.stats import spearmanr
import xlsxwriter

from fasta import read_records, write_records, Record
from utilities import async_call, safe_phylo_read, safe_delete, _self_path, safe_mkdir, chunks, make_bold_formatting, \
    make_p_formatting, make_rho_formatting

ErcResult = namedtuple('ErcResult', ['file1', 'file2', 'rho', 'p', 'raw_rates'])

ErcDataSource = namedtuple('ErcDataSource', ['tree_dir', 'ercs', 'sep', 'is_oneway'])
ErcDataSourceEntry = namedtuple('ErcDataSourceEntry', ['path1', 'path2', 'rho', 'p'])
SegmentedERC = namedtuple("SegmentedERC", ["result", 'is_file1_segmented', 'is_file2_segmented',
                                           'file1_base', 'file2_base', 'file1_range', 'file2_range'])

datasources = []
GAP_CHARS = {'?', '-', '*'}

_10mya_cutoff = """ORNITHORHYNCHUS_ANATINUS
MONODELPHIS_DOMESTICA
PHASCOLARCTOS_CINEREUS
SARCOPHILUS_HARRISII
DASYPUS_NOVEMCINCTUS
LOXODONTA_AFRICANA
TRICHECHUS_MANATUS
ORYCTEROPUS_AFER
ELEPHANTULUS_EDWARDII
CHRYSOCHLORIS_ASIATICA
ECHINOPS_TELFAIRI
OCHOTONA_PRINCEPS
ORYCTOLAGUS_CUNICULUS
FUKOMYS_DAMARENSIS
HETEROCEPHALUS_GLABER
CAVIA_PORCELLUS
CHINCHILLA_LANIGERA
OCTODON_DEGUS
MARMOTA_MARMOTA
CASTOR_CANADENSIS
DIPODOMYS_ORDII
JACULUS_JACULUS
NANNOSPALAX_GALILI
PEROMYSCUS_MANICULATUS
MICROTUS_OCHROGASTER
CRICETULUS_GRISEUS
MESOCRICETUS_AURATUS
MERIONES_UNGUICULATUS
RATTUS_NORVEGICUS
MUS_MUSCULUS
GALEOPTERUS_VARIEGATUS
OTOLEMUR_GARNETTII
MICROCEBUS_MURINUS
PROPITHECUS_COQUERELI
CARLITO_SYRICHTA
AOTUS_NANCYMAAE
CALLITHRIX_JACCHUS
SAIMIRI_BOLIVIENSIS
NOMASCUS_LEUCOGENYS
PONGO_ABELII
HOMO_SAPIENS
GORILLA_GORILLA
RHINOPITHECUS_BIETI
RHINOLOPHUS_ROXELLANA
COLOBUS_ANGOLENSIS
PILIOCOLOBUS_TEPHROSCELES
CHLOROCEBUS_SABAEUS
PAPIO_ANUBIS
MANDRILLUS_LEUCOPHAEUS
MACACA_MULATTA
CONDYLURA_CRISTATA
ERINACEUS_EUROPAEUS
SOREX_ARANEUS
HIPPOSIDEROS_ARMIGER
RHINOLOPHUS_SINICUS
ROUSETTUS_AEGYPTIACUS
PTEROPUS_ALECTO
PTEROPUS_VAMPYRUS
MINIOPTERUS_NATALENSIS
EPTESICUS_FUSCUS
MYOTIS_DAVIDII
MYOTIS_BRANDTII
MYOTIS_LUCIFUGUS
CERATOTHERIUM_SIMUM
EQUUS_CABALLUS
MANIS_JAVANICA
ACINONYX_JUBATUS
FELIS_CATUS
CANIS_LUPUS
AILUROPODA_MELANOLEUCA
URSUS_MARITIMUS
ENHYDRA_LUTRIS
MUSTELA_PUTORIUS
LEPTONYCHOTES_WEDDELLII
ODOBENUS_ROSMARUS
VICUGNA_PACOS
CAMELUS_FERUS
SUS_SCROFA
BALAENOPTERA_ACUTOROSTRATA
PHYSETER_CATODON
LIPOTES_VEXILLIFER
DELPHINAPTERUS_LEUCAS
ORCINUS_ORCA
ODOCOILEUS_VIRGINIANUS
PANTHOLOPS_HODGSONII
CAPRA_HIRCUS
OVIS_ARIES
BUBALUS_BUBALIS
BOS_GRUNNIENS""".splitlines()

_20mya_cutoff = """ORNITHORHYNCHUS_ANATINUS
MONODELPHIS_DOMESTICA
PHASCOLARCTOS_CINEREUS
SARCOPHILUS_HARRISII
DASYPUS_NOVEMCINCTUS
LOXODONTA_AFRICANA
TRICHECHUS_MANATUS
ORYCTEROPUS_AFER
ELEPHANTULUS_EDWARDII
CHRYSOCHLORIS_ASIATICA
ECHINOPS_TELFAIRI
OCHOTONA_PRINCEPS
ORYCTOLAGUS_CUNICULUS
FUKOMYS_DAMARENSIS
HETEROCEPHALUS_GLABER
CAVIA_PORCELLUS
CHINCHILLA_LANIGERA
OCTODON_DEGUS
MARMOTA_MARMOTA
CASTOR_CANADENSIS
DIPODOMYS_ORDII
JACULUS_JACULUS
NANNOSPALAX_GALILI
PEROMYSCUS_MANICULATUS
MICROTUS_OCHROGASTER
CRICETULUS_GRISEUS
MESOCRICETUS_AURATUS
MERIONES_UNGUICULATUS
RATTUS_NORVEGICUS
MUS_MUSCULUS
GALEOPTERUS_VARIEGATUS
OTOLEMUR_GARNETTII
MICROCEBUS_MURINUS
PROPITHECUS_COQUERELI
CARLITO_SYRICHTA
AOTUS_NANCYMAAE
CALLITHRIX_JACCHUS
SAIMIRI_BOLIVIENSIS
NOMASCUS_LEUCOGENYS
HOMO_SAPIENS
RHINOPITHECUS_BIETI
MACACA_MULATTA
CONDYLURA_CRISTATA
ERINACEUS_EUROPAEUS
SOREX_ARANEUS
HIPPOSIDEROS_ARMIGER
RHINOLOPHUS_SINICUS
ROUSETTUS_AEGYPTIACUS
PTEROPUS_VAMPYRUS
MINIOPTERUS_NATALENSIS
EPTESICUS_FUSCUS
MYOTIS_DAVIDII
MYOTIS_LUCIFUGUS
CERATOTHERIUM_SIMUM
EQUUS_CABALLUS
MANIS_JAVANICA
FELIS_CATUS
CANIS_LUPUS
AILUROPODA_MELANOLEUCA
URSUS_MARITIMUS
ENHYDRA_LUTRIS
MUSTELA_PUTORIUS
LEPTONYCHOTES_WEDDELLII
ODOBENUS_ROSMARUS
VICUGNA_PACOS
CAMELUS_FERUS
SUS_SCROFA
BALAENOPTERA_ACUTOROSTRATA
PHYSETER_CATODON
LIPOTES_VEXILLIFER
DELPHINAPTERUS_LEUCAS
ORCINUS_ORCA
ODOCOILEUS_VIRGINIANUS
PANTHOLOPS_HODGSONII
OVIS_ARIES
BOS_GRUNNIENS""".splitlines()

_30mya_cutoff = """ORNITHORHYNCHUS_ANATINUS
MONODELPHIS_DOMESTICA
PHASCOLARCTOS_CINEREUS
SARCOPHILUS_HARRISII
DASYPUS_NOVEMCINCTUS
LOXODONTA_AFRICANA
TRICHECHUS_MANATUS
ORYCTEROPUS_AFER
ELEPHANTULUS_EDWARDII
CHRYSOCHLORIS_ASIATICA
ECHINOPS_TELFAIRI
OCHOTONA_PRINCEPS
ORYCTOLAGUS_CUNICULUS
FUKOMYS_DAMARENSIS
HETEROCEPHALUS_GLABER
CAVIA_PORCELLUS
CHINCHILLA_LANIGERA
OCTODON_DEGUS
MARMOTA_MARMOTA
CASTOR_CANADENSIS
DIPODOMYS_ORDII
JACULUS_JACULUS
NANNOSPALAX_GALILI
PEROMYSCUS_MANICULATUS
MESOCRICETUS_AURATUS
MERIONES_UNGUICULATUS
MUS_MUSCULUS
GALEOPTERUS_VARIEGATUS
OTOLEMUR_GARNETTII
MICROCEBUS_MURINUS
PROPITHECUS_COQUERELI
CARLITO_SYRICHTA
CALLITHRIX_JACCHUS
HOMO_SAPIENS
MACACA_MULATTA
CONDYLURA_CRISTATA
ERINACEUS_EUROPAEUS
SOREX_ARANEUS
HIPPOSIDEROS_ARMIGER
RHINOLOPHUS_SINICUS
ROUSETTUS_AEGYPTIACUS
PTEROPUS_VAMPYRUS
MINIOPTERUS_NATALENSIS
EPTESICUS_FUSCUS
MYOTIS_LUCIFUGUS
CERATOTHERIUM_SIMUM
EQUUS_CABALLUS
MANIS_JAVANICA
FELIS_CATUS
CANIS_LUPUS
AILUROPODA_MELANOLEUCA
MUSTELA_PUTORIUS
ODOBENUS_ROSMARUS
CAMELUS_FERUS
SUS_SCROFA
BALAENOPTERA_ACUTOROSTRATA
PHYSETER_CATODON
ORCINUS_ORCA
ODOCOILEUS_VIRGINIANUS
BOS_GRUNNIENS""".splitlines()


def read_gene_info_as_l2n(filename: str, as_symbols: bool = False) -> Dict[str, str]:
    """
    Reads a gene info list file to make a label -> name dict.
    :param filename: The file.
    :param as_symbols: Whether to return symbols or full names. If none, both are combined
    :return: The l2n dict.
    """
    l2n = {
        "209332at40674": "APOA1" if as_symbols else ("Apolipoprotein A-I (APOA1)" if as_symbols is None else "Apolipoprotein A-I"),
        "91651at40674": "CETP" if as_symbols else ("Cholesteryl ester transfer protein (CETP)" if as_symbols is None else "Cholesteryl ester transfer protein"),
        "209652at40674": "APOC1" if as_symbols else ("Apolipoprotein C-I (APOC1)" if as_symbols is None else "Apolipoprotein C-I"),
        "158878at40674": "VEGFA" if as_symbols else ("Vascular endothelial growth factor A (VEGFA)" if as_symbols is None else "Vascular endothelial growth factor A"),
        "coag9": "F9" if as_symbols else ("Coagulation Factor 9 (F9)" if as_symbols is None else "Coagulation Factor 9"),
        "coag10": "F10" if as_symbols else ("Coagulation Factor 10 (F10)" if as_symbols is None else "Coagulation Factor 10"),
        "68161at40674": "POLH" if as_symbols else ("DNA polymerase eta (POLH)" if as_symbols is None else "DNA polymerase eta"),
        "71251at40674": "SELE" if as_symbols else ("Selectin E (SELE)" if as_symbols is None else "Selectin E"),
        "147514at40674": "ALKBH4" if as_symbols else ("alkB homolog 4, lysine demethylase (ALKBH4)" if as_symbols is None else "alkB homolog 4, lysine demethylase"),
        "70809at40674": "MARCHF10" if as_symbols else ("Testis secretory sperm-binding protein Li 228n (ALKBH4)" if as_symbols is None else "Testis secretory sperm-binding protein Li 228n")
    }  # FIXME: Update the master file to include these builtins
    with open(filename, 'r') as f:
        first = True
        for l in f:
            if first:
                first = False
                continue
            row = l.strip().split("\t")
            if len(row) == 0:
                continue
            offset = 1 if as_symbols else 0
            odb = row[0]

            if odb in l2n:
                continue

            human = row[3 + offset]
            mouse = row[5 + offset]
            chimp = row[7 + offset]
            bonobo = row[9 + offset]
            generic = row[11]
            entrez = row[12]
            uniprot = row[14]

            name = 'NOT FOUND'

            if human == 'NOT FOUND':
                if mouse == 'NOT FOUND':
                    if chimp == 'NOT FOUND':
                        if bonobo == 'NOT FOUND':
                            if entrez == 'NOT FOUND' or not as_symbols:
                                if uniprot == 'NOT FOUND' or as_symbols:
                                    name = odb if as_symbols else generic
                                else:
                                    name = uniprot
                            else:
                                name = entrez
                        else:
                            name = bonobo
                    else:
                        name = chimp
                else:
                    name = mouse
            else:
                name = human

            l2n[odb] = name

            if as_symbols is None:
                symbol = "NOT FOUND"
                if row[4] == 'NOT FOUND':
                    if row[6] == 'NOT FOUND':
                        if row[8] == 'NOT FOUND':
                            symbol = row[10]
                        else:
                            symbol = row[8]
                    else:
                        symbol = row[6]
                else:
                    symbol = row[4]
                l2n[odb] = f"{l2n[odb]} ({symbol})"

    return l2n


def make_l2n(as_symbols = None) -> Dict[str, str]:
    annotations = osp.join(_self_path(), "data", "full_gene_info.tsv") # Latest matrix
    return read_gene_info_as_l2n(annotations, as_symbols)


async def get_rates(timetree: Union[str, PhyloTree], prune: bool = True, taxa: List[str] = None, *trees: Union[str, PhyloTree]) -> Tuple[List[str], List[List[float]]]:
    """
    Get rates from time + trees.
    :param timetree: The time tree.
    :param prune: Whether the prune the tree based on the passed in taxa. Trees are still pruned to match species.
    :param taxa: The taxa to consider.
    :param trees: The trees to get rates for.
    :return: A tuple, first you get a list representing the taxon indices. Second, you get a list of lists. The outer
        list represents the tree (note: index 0 = time), the inner list corresponds to taxon rates.
    """
    timetree = timetree if not isinstance(timetree, str) else safe_phylo_read(timetree)
    trees = [tree if not isinstance(tree, str) else safe_phylo_read(tree) for tree in trees]

    shared_taxa = set(timetree.iter_leaf_names())
    for tree in trees:
        shared_taxa = shared_taxa & set(tree.iter_leaf_names())

    if taxa:
        taxa = list(shared_taxa & set(taxa))

    shared_taxa = list(shared_taxa)

    if not taxa:
        taxa = shared_taxa

    if prune:
        timetree = prune_tree(timetree, taxa)
        trees = [(prune_tree(tree, taxa)) for tree in trees]
    else:
        timetree = prune_tree(timetree, shared_taxa)
        trees = [(prune_tree(tree, shared_taxa)) for tree in trees]

    rates = [[] for _ in trees]
    rates.append([])
    leaves = []
    for leaf in timetree.iter_leaves():
        if leaf.name not in taxa:
            continue

        leaves.append(leaf.name)
        time = np.float(leaf.get_distance(leaf.up))
        rates[0].append(time)
        for i, tree in enumerate(trees):
            if time == 0.:
                rates[i+1].append(np.float(0.0))
            else:
                leaf2 = list(tree.get_leaves_by_name(leaf.name))[0]
                rates[i+1].append(np.float(leaf2.get_distance(leaf2.up)) / time)

    return leaves, rates


async def get_changed_branches(timetree: PhyloTree, tree: PhyloTree, prune_list: List[str]) -> List[str]:
    pruned_time_tree = prune_tree(timetree, list(tree.iter_leaf_names()))
    branches = {leaf.name: leaf.get_distance(leaf.up) for leaf in pruned_time_tree.iter_leaves()}
    pruned = prune_tree(pruned_time_tree, prune_list)
    changed = []
    for leaf in pruned.iter_leaves():
        if branches[leaf.name] < leaf.get_distance(leaf.up):
            changed.append(leaf.name)
    return changed


# Convert rate object to correlation
def rates_to_correlation(res: Tuple[List[str], List[List[float]]], corr_to_time: bool = False) -> Tuple[np.float, np.float]:
    results = spearmanr(np.array(res[1][0 + (not corr_to_time)]), np.array(res[1][1 + (not corr_to_time)]), nan_policy='raise')

    return results.correlation, results.pvalue


async def archive_directory(directory: str, archive: str):
    """
    Archives a directory.

    :param directory: The directory to archive
    :param archive: The archive name
    """
    if osp.exists(directory) and not osp.exists(archive):
        try:
            await async_call(f"tar -cjf {archive} {directory}")
            shutil.rmtree(directory, True)
        except Exception as e:
            print(f"WARNING: Unable to tar compress the directory {directory}")
            traceback.print_exc()


async def clean_fasta(treefile: str, fastafile: str, clean_out: str = None):
    """
    Given a topology, normalizes fasta taxa naming and drops all taxa with paralogs.

    :param treefile: The tree topology.
    :param fastafile: The input file.
    :param clean_out: The output file path.
    """
    if not clean_out:
        clean_out = fastafile

    records = read_records(fastafile)

    taxa_names = set(r.title for r in records)
    fasta2tree = dict()
    for n in safe_phylo_read(treefile).get_leaf_names():
        norm_n = n.replace(' ', '_').upper().strip()
        for t in taxa_names:
            norm_t = t.replace(' ', '_').upper().strip()
            if norm_n in norm_t or norm_t in norm_n:
                fasta2tree[t] = n

    counts = dict()
    for r in records:
        r.title = fasta2tree.get(r.title, r.title)
        if r.title not in counts:
            counts[r.title] = 0
        counts[r.title] += 1
    dupes = {taxon for taxon, count in counts.items() if count > 1}

    new_records = []
    for r in records:
        if r.title not in dupes:
            new_records.append(r)

    write_records(clean_out, new_records)


async def align(fastafile: str, outputfile: str):
    """
    Aligns a file.

    :param fastafile: The input path
    :param outputfile: The output path
    """
    assert fastafile and osp.exists(fastafile)

    if not outputfile:
        outputfile = fastafile

    # NOTE: --anysymbol allows the U symbol to be used in protein seqs https://mafft.cbrc.jp/alignment/software/anysymbol.html
    await async_call(f"mafft --maxiterate 1000 --localpair --anysymbol --thread 5 --out {outputfile} {fastafile}")

    assert osp.exists(outputfile)


async def trim(inputfile: str, outputfile: str):
    """
    Trims an alignment.

    :param inputfile: The input alignment.
    :param outputfile: The trimmed alignment.
    :return:
    """
    assert inputfile and osp.exists(inputfile)

    await async_call("trimal -in {} -out {} -automated1".format(inputfile, outputfile))


def prune_tree(in_file: Union[str, PhyloTree], retained: List[str], out_group: Optional[str] = None,
               out_file: str = None, deroot: bool = False) -> Optional[PhyloTree]:

    if isinstance(in_file, str):
        assert in_file and osp.exists(in_file)

        if out_file and osp.exists(out_file):
            os.remove(out_file)

        tree = safe_phylo_read(in_file)
    else:
        tree = in_file.copy()

    curr_leaves = [l for l in tree.iter_leaf_names()]

    try:
        tree.prune([r for r in retained if r in curr_leaves], preserve_branch_length=True)
    except Exception as e:
        print(in_file)
        print(retained)
        print(curr_leaves)
        print([r for r in retained if r in curr_leaves])
        raise e
    tree.resolve_polytomy()

    if out_group is not None:
        tree.set_outgroup(out_group)

    if deroot:
        tree.unroot("keep")

    if out_file:
        tree.write(outfile=out_file)
    else:
        return tree


async def prune(alignment: str, tree: str, treeout: str):
    """
    Prunes a tree based on taxa in an alignment.

    :param alignment: The alignment.
    :param tree: The tree.
    :param treeout: The output path.
    """
    prune_tree(tree, [r.title for r in read_records(alignment) if '.' not in r.title], out_group=None, out_file=treeout)


async def iqtree2(alignment: str,
                  consensus_tree: str = None,
                  model: str = "AUTO",  # LG+F+G+I
                  output: str = None,
                  bootstraps: int = 1000,
                  cpus: int = 2,
                  expect_model_violations: bool = False,
                  partitions_file: str = None,
                  site_specific_rates: bool = True,
                  clean_outputs: bool = True,
                  ancestral_states: bool = True):
    # Outputs: {alignment}.iqtree -> full report
    #   {alignment}.treefile -> tree
    #   {alignment}.log -> log file
    #   {alignment}.model -> selected model if none is specified
    #   {alignment}.contree -> tree with bootstrap values only
    #   {alignment}.splits.nex -> support values in percentage for all splits (bipartitions), computed as the occurence frequencies in the bootstrap trees.
    #   {alignment}.rate -> site specific rates. See http://www.iqtree.org/doc/Advanced-Tutorial

    cmd = f"iqtree2 -s {alignment} -B {bootstraps} -T {cpus} -st AA -seed 1234567890"

    if model != "AUTO":
        cmd += f" -m {model}"  # GTR+I+G

    if expect_model_violations:  # . With this option UFBoot will further optimize each bootstrap tree using a hill-climbing nearest neighbor interchange (NNI) search based directly on the corresponding bootstrap alignment.
        cmd += " -bnni"

    if partitions_file:
        cmd += f" -p {partitions_file}"

    if consensus_tree:
        cmd += f" -g {consensus_tree}"

    if site_specific_rates:
        cmd += " --mlrate"

    if ancestral_states:
        cmd += " -asr"

    await async_call(cmd)

    if output:
        base = partitions_file if partitions_file else alignment
        shutil.move(base + ".treefile", output)

        if site_specific_rates:
            shutil.move(base + ".mlrate", output + ".mlrate")

        if ancestral_states:
            shutil.move(base + ".state", output + ".state")

    if clean_outputs:
        safe_delete(base + ".ckp.gz")  # Checkpoint info
        safe_delete(base + ".log")  # Log
        safe_delete(base + ".uniqueseq.phy")  # Pruned alignment
        safe_delete(base + ".contree")  # consensus tree
        safe_delete(base + ".iqtree")  # more logging?
        safe_delete(base + ".parstree")  # Parsimony tree
        safe_delete(base + ".splits.nex")  # Split info


async def protein_tree(alignment: str, topology: str, output: str, partition_file: str = None):
    """
    Generates a tree using iqtree.

    :param alignment: The alignment.
    :param topology: The topology.
    :param output: The output file. Site specific rates: output + ".mlrate" and ancestral states: output + ".state"
    """
    try:
        if not output or not osp.exists(output):
            await iqtree2(alignment, topology, model="LG+F+G+I", output=output, ancestral_states=True, partitions_file=partition_file)
    except Exception as e:
        print(f"WARNING: Could not generate tree for {alignment}!", flush=True)
        traceback.print_exc()


async def partial_correlation(tree1: PhyloTree, tree2: PhyloTree, to_control: PhyloTree, species: PhyloTree,
                        method: str = "spearman", correct_control: bool = True,
                        to_control2: PhyloTree = None, correct_control2: bool = True,
                        to_control3: PhyloTree = None, correct_control3: bool = True,
                        to_control4: PhyloTree = None, correct_control4: bool = True,
                        normalize_controls: bool = False, return_shared_taxa_count: bool = False) -> Union[Tuple[np.float, np.float], Tuple[np.float, np.float, np.int]]:
    """
    Bound using https://cran.r-project.org/web/packages/ppcor/ppcor.pdf

    Requires R

    returns corr, p-val
    """

    fingerprint = str(hex((hash(tree1) ^ hash(tree2) ^ hash(species)) + hash(to_control) + hash(method) \
                          + hash(to_control2) + hash(to_control3) + hash(to_control4) + hash(correct_control) + hash(correct_control2) \
                          + hash(correct_control3) + hash(correct_control4) + hash(normalize_controls)))

    if not to_control2:
        all_taxa = list(set(tree1.iter_leaf_names()) & set(tree2.iter_leaf_names()) & set(to_control.iter_leaf_names()))
    elif not to_control3:
        all_taxa = list(set(tree1.iter_leaf_names()) & set(tree2.iter_leaf_names()) & set(to_control2.iter_leaf_names()) & set(to_control.iter_leaf_names()))
    elif not to_control4:
        all_taxa = list(set(tree1.iter_leaf_names()) & set(tree2.iter_leaf_names()) & set(to_control3.iter_leaf_names()) & set(to_control2.iter_leaf_names()) & set(to_control.iter_leaf_names()))
    else:
        all_taxa = list(set(tree1.iter_leaf_names()) & set(tree2.iter_leaf_names()) & set(to_control4.iter_leaf_names()) & set(to_control3.iter_leaf_names()) & set(to_control2.iter_leaf_names()) & set(to_control.iter_leaf_names()))

    tree1 = tree1.copy()
    tree2 = tree2.copy()
    to_control = to_control.copy()
    to_control2 = to_control2.copy() if to_control2 else None
    to_control3 = to_control3.copy() if to_control3 else None
    to_control4 = to_control4.copy() if to_control4 else None

    species = species.copy()

    tree1.prune(all_taxa, preserve_branch_length=True)
    tree2.prune(all_taxa, preserve_branch_length=True)
    to_control.prune(all_taxa, preserve_branch_length=True)
    species.prune(all_taxa, preserve_branch_length=True)
    if to_control2:
        to_control2.prune(all_taxa, preserve_branch_length=True)
    if to_control3:
        to_control3.prune(all_taxa, preserve_branch_length=True)
    if to_control4:
        to_control4.prune(all_taxa, preserve_branch_length=True)

    tree1_data = dict()
    tree2_data = dict()
    control_data = dict()
    species_data = dict()
    control_data2 = dict()
    control_data3 = dict()
    control_data4 = dict()

    for leaf in tree1.iter_leaves():
        if leaf.name not in all_taxa:
            continue
        tree1_data[leaf.name] = leaf.get_distance(leaf.up)
    for leaf in tree2.iter_leaves():
        if leaf.name not in all_taxa:
            continue
        tree2_data[leaf.name] = leaf.get_distance(leaf.up)
    for leaf in to_control.iter_leaves():
        if leaf.name not in all_taxa:
            continue
        control_data[leaf.name] = leaf.get_distance(leaf.up)
    for leaf in species.iter_leaves():
        if leaf.name not in all_taxa:
            continue
        species_data[leaf.name] = leaf.get_distance(leaf.up)
    if to_control2:
        for leaf in to_control2.iter_leaves():
            if leaf.name not in all_taxa:
                continue
            control_data2[leaf.name] = leaf.get_distance(leaf.up)
    if to_control3:
        for leaf in to_control3.iter_leaves():
            if leaf.name not in all_taxa:
                continue
            control_data3[leaf.name] = leaf.get_distance(leaf.up)
    if to_control4:
        for leaf in to_control4.iter_leaves():
            if leaf.name not in all_taxa:
                continue
            control_data4[leaf.name] = leaf.get_distance(leaf.up)

    if normalize_controls:
        max_control1 = max(control_data.values())
        min_control1 = min(control_data.values())
        for k, v in list(control_data.items()):
            control_data[k] = (v - min_control1) / (max_control1 - min_control1)

        if to_control2:
            max_control2 = max(control_data2.values())
            min_control2 = min(control_data2.values())
            for k, v in list(control_data2.items()):
                control_data2[k] = (v - min_control2) / (max_control2 - min_control2)
        if to_control3:
            max_control3 = max(control_data3.values())
            min_control3 = min(control_data3.values())
            for k, v in list(control_data3.items()):
                control_data3[k] = (v - min_control3) / (max_control3 - min_control3)
        if to_control4:
            max_control4 = max(control_data4.values())
            min_control4 = min(control_data4.values())
            for k, v in list(control_data4.items()):
                control_data4[k] = (v - min_control4) / (max_control4 - min_control4)

    with open(f"temp_{fingerprint}.csv", 'w') as f:
        f.write("x,y,")
        if to_control4:
            f.write("z1,z2,z3,z4\n")
        elif to_control3:
            f.write("z1,z2,z3\n")
        elif to_control2:
            f.write("z1,z2\n")
        else:
            f.write("z\n")
        for taxon in all_taxa:
            f.write(f"{np.float(tree1_data[taxon]) / np.float(species_data[taxon])},")
            f.write(f"{np.float(tree2_data[taxon]) / np.float(species_data[taxon])},")
            if correct_control:
                f.write(f"{np.float(control_data[taxon]) / np.float(species_data[taxon])}")
            else:
                f.write(f"{np.float(control_data[taxon])}")

            if to_control2:
                if correct_control2:
                    f.write(f",{np.float(control_data2[taxon]) / np.float(species_data[taxon])}")
                else:
                    f.write(f",{np.float(control_data2[taxon])}")

            if to_control3:
                if correct_control3:
                    f.write(f",{np.float(control_data3[taxon]) / np.float(species_data[taxon])}")
                else:
                    f.write(f",{np.float(control_data3[taxon])}")

            if to_control4:
                if correct_control4:
                    f.write(f",{np.float(control_data4[taxon]) / np.float(species_data[taxon])}")
                else:
                    f.write(f",{np.float(control_data4[taxon])}")

            f.write("\n")

    if to_control4:
        script_name = "specific_partial_correlation4.R"
    elif to_control3:
        script_name = "specific_partial_correlation3.R"
    elif to_control2:
        script_name = "specific_partial_correlation2.R"
    else:
        script_name = "specific_partial_correlation.R"

    await async_call(f"Rscript --vanilla {osp.join(_self_path(), 'R', script_name)} temp_{fingerprint}.csv {method} temp_out_{fingerprint}.csv")

    nums = []
    if not osp.exists(f'temp_out_{fingerprint}.csv'):
        nums.append(np.nan)
        nums.append(np.nan)
    else:
        with open(f'temp_out_{fingerprint}.csv', 'r') as f:
            first = True
            for l in f:
                if first:
                    first = False
                    continue
                if len(l.strip()) == 0:
                    continue
                num_str = l.strip()
                if num_str == 'NA':
                    nums.append(np.nan)
                else:
                    nums.append(np.float(num_str))

    os.remove(f"temp_{fingerprint}.csv")
    os.remove(f"temp_out_{fingerprint}.csv")

    if return_shared_taxa_count:
        return nums[0], nums[1], np.int(len(all_taxa))
    return nums[0], nums[1]


async def correlate(tree1: Union[str, PhyloTree], tree2: Union[str, PhyloTree], norm_tree: Union[str, PhyloTree],
                    internal_requirement: float = -1, include_terminal: bool = False,
                    ignore_species: List[str] = None) -> Tuple[float, float, Dict[str, Tuple[float, float, float]]]:
    norm_tree = safe_phylo_read(norm_tree) if isinstance(norm_tree, str) else norm_tree.copy()  # We can reasonably assume this won't error
    try:
        tree1 = safe_phylo_read(tree1) if isinstance(tree1, str) else tree1.copy()
        tree2 = safe_phylo_read(tree2) if isinstance(tree2, str) else tree2.copy()
    except Exception as e:
        print(f"Problem handling trees: {tree1}, {tree2}, {norm_tree}", file=sys.stderr)
        return nan, nan, {l: (nan, nan) for l in norm_tree.iter_leaves()}

    tree1_rates = list()
    tree2_rates = list()

    raw_rates = {}

    if internal_requirement < 0:
        for leaf in tree1.iter_leaves():
            if ignore_species and leaf.name in ignore_species:
                continue

            tree1_rate = leaf.get_distance(leaf.up)
            results = tree2.get_leaves_by_name(leaf.name)
            if leaf.name not in [x.name for x in results]:
                print(f"Skipping {leaf.name}...", flush=True)
                continue
            tree2_leaf = results[0]
            assert tree2_leaf.is_leaf()
            tree2_rate = tree2_leaf.get_distance(tree2_leaf.up)
            results = norm_tree.get_leaves_by_name(leaf.name)
            if leaf.name not in [x.name for x in results]:
                print(f"Skipping {leaf.name}...", flush=True)
                continue
            norm_leaf = results[0]
            assert norm_leaf.is_leaf()
            norm_factor = norm_leaf.get_distance(norm_leaf.up)
            tree1_rates.append(0 if norm_factor == 0. else tree1_rate / norm_factor)
            tree2_rates.append(0 if norm_factor == 0. else tree2_rate / norm_factor)
            raw_rates[leaf.name] = (tree1_rates[-1], tree2_rates[-1], norm_factor)
            # else:
            #     raise AssertionError("Unexpected norm_mode " + norm_mode)
    else:
        for node in tree1.traverse():
            if node.is_leaf():
                if not include_terminal or (ignore_species and node.name in ignore_species):
                    continue
                leaf = node
                tree1_rate = leaf.get_distance(leaf.up)
                results = tree2.get_leaves_by_name(leaf.name)
                if leaf.name not in [x.name for x in results]:
                    print(f"Skipping {leaf.name}...", flush=True)
                    continue
                tree2_leaf = results[0]
                assert tree2_leaf.is_leaf()
                tree2_rate = tree2_leaf.get_distance(tree2_leaf.up)
                results = norm_tree.get_leaves_by_name(leaf.name)
                if leaf.name not in [x.name for x in results]:
                    print(f"Skipping {leaf.name}...", flush=True)
                    continue
                norm_leaf = results[0]
                assert norm_leaf.is_leaf()
                norm_factor = norm_leaf.get_distance(norm_leaf.up)
                tree1_rates.append(0 if norm_factor == 0. else tree1_rate / norm_factor)
                tree2_rates.append(0 if norm_factor == 0. else tree2_rate / norm_factor)
                raw_rates[leaf.name] = (tree1_rates[-1], tree2_rates[-1], norm_factor)
                # else:
                #     raise AssertionError("Unexpected norm_mode " + norm_mode)
            else:
                leaves = set(node.get_leaf_names())

                tree2_node = None
                for node2 in tree2.traverse():
                    if leaves == set(node2.get_leaf_names()):
                        tree2_node = node2
                        break

                if not tree2_node:
                    print(f"Skipping {node.name}...", flush=True)
                    continue

                norm_node = None
                for nnode in norm_tree.traverse():
                    if leaves == set(nnode.get_leaf_names()):
                        norm_node = nnode
                        break

                if not norm_node:
                    print(f"Skipping {node.name}...", flush=True)
                    continue

                if norm_node.up is not None:
                    if norm_node.get_distance(norm_node.up) < internal_requirement:
                        continue

                    norm_factor = norm_node.get_distance(norm_node.up)

                    tree1_rate = node.get_distance(node.up)
                    tree2_rate = tree2_node.get_distance(tree2_node.up)

                    tree1_rates.append(0 if norm_factor == 0. else tree1_rate / norm_factor)
                    tree2_rates.append(0 if norm_factor == 0. else tree2_rate / norm_factor)
                    raw_rates["$".join(leaves)] = (tree1_rates[-1], tree2_rates[-1], norm_factor)

    results = spearmanr(array(tree1_rates), array(tree2_rates), nan_policy='raise')

    return results.correlation, results.pvalue, raw_rates


async def generate_erc(topology: str, alignment1: str, alignment2: str, internal_requirement: float = -1,
                       include_terminal: bool = False, taxon_set: List[str] = None,
                       time_corrected: bool = False) -> ErcResult:
    """
    Calculates an ERC.

    :param topology: The fixed topology file.
    :param alignment1: The tree for alignment 1.
    :param alignment2: The tree for alignment 2.
    :param internal_requirement: Threshold for correlating internal branches (-1 = no correlation)
    :param include_terminal: Whether to include the terminal branches with internal ERC calculations
    :param taxon_set: The taxa to consider.
    :param time_corrected: Whether to run time corrected pcorrs.
    :return: The result.
    """

    if not osp.exists(alignment1) or not osp.exists(alignment2):
        return ErcResult(alignment1, alignment2, 0, 1, dict())

    shared_taxa = {l for l in safe_phylo_read(alignment1).iter_leaf_names()} \
                  & {l for l in safe_phylo_read(alignment2).iter_leaf_names()}

    if taxon_set:
        shared_taxa = shared_taxa & set(taxon_set)

    shared_taxa = list(shared_taxa)

    prune_top, prune_aln1, prune_aln2 = await asyncio.gather(prune_tree(topology, shared_taxa, None),
                                                             prune_tree(alignment1, shared_taxa, None),
                                                             prune_tree(alignment2, shared_taxa, None))

    if time_corrected:
        rho, p = await partial_correlation(prune_aln1, prune_aln2, prune_top, prune_top,
                                           correct_control=False)
        raw = None
    else:
        rho, p, raw = await correlate(prune_aln1, prune_aln2, prune_top,
                                      internal_requirement=internal_requirement, include_terminal=include_terminal)

    return ErcResult(alignment1, alignment2, rho, p, raw)


def register_erc_datasource(treedir: str, datafile: str, sep: str = ',', is_oneway_to: str = None):
    """
    Registers a data source for pre-genned ercs and trees.

    :param treedir: A directory containing trees referred to in the data file.
    :param datafile: The data file with columns for alignment1, alignment2, erc, p
    :param sep: The field separator.
    :param is_oneway_to: If there is only one file column, what is the implicit file.
    """
    datasources.append(ErcDataSource(treedir, datafile, sep, is_oneway_to))


def clear_datasources():
    """
    Clears the registered data sources.
    """
    global datasources
    datasources.clear()


def read_datasource(datasource: ErcDataSource) -> List[ErcDataSourceEntry]:
    """
    Reads a data source.

    :param datasource: The data source to read.
    :return: The entries in the data source.
    """
    entries = []
    first = True
    with open(datasource.ercs, 'r') as f:
        for l in f:
            if first:
                first = False
                continue
            try:
                split = l.strip().split(datasource.sep)
                align1 = split[0]
                if datasource.is_oneway:
                    align2 = datasource.is_oneway
                    offset = 0
                else:
                    align2 = split[1]
                    offset = 1
                erc = np.float(split[1 + offset])
                p = np.float(split[2 + offset])

                entries.append(ErcDataSourceEntry(osp.join(datasource.tree_dir, align1),
                                                  osp.join(datasource.tree_dir, align2),
                                                  erc, p))
            except:
                pass

    return entries


async def generate_network() -> nx.Graph:
    """
    Builds a network from currently registered data sources.
    :return: The erc network.
    """
    graph = nx.Graph()
    for datasource in datasources:
        for entry in read_datasource(datasource):
            id1 = osp.basename(entry.path1).split('.')[0]
            id2 = osp.basename(entry.path2).split('.')[0]

            graph.add_edge(id1, id2, rho=entry.rho, p=entry.p, weight=entry.rho)
    return graph


def enrich_network(net: List[str], background_net: List[str], id2Name: Dict[str, str], basefile: str):
    """
    Runs an enrichment on a network of genes.
    :param net: The network nodes to check for enrichment.
    :param background_net: The background network nodes for the enrichment.
    :param id2Name: A mapping from node id -> unambiguous gene symbol.
    :param basefile: The base name for results.
    """
    background = []
    gene_list = []
    for n in net:
        if id2Name.get(n, n) != 'NOT FOUND':
            for split in id2Name.get(n, n).split(";"):
                gene_list.append(split.split("_")[-1])
    for n in background_net:
        if id2Name.get(n, n) != 'NOT FOUND':
            for split in id2Name.get(n, n).split(";"):
                background.append(split.split("_")[-1])

    if len([g for g in gene_list if 'NOT FOUND' not in g]) < 2:
        return

    with open(basefile + ".gene_list.txt", 'w') as f:
        for g in gene_list:
            f.write(g)
            f.write("\n")

    kegg_enr = gp.enrichr(gene_list,
                     'KEGG_2019_Human',
                     'human',
                     'kegg_output',
                     osp.join(basefile, 'kegg_output'),
                     background,
                     no_plot=True)

    kegg_enr.results.to_csv(basefile + '.kegg.csv')
    shutil.rmtree(osp.join(basefile))

    go_enr = gp.enrichr(gene_list,
                          'GO_Biological_Process_2018',
                          'human',
                          'go_output1',
                          osp.join(basefile, 'go_output'),
                          background,
                          no_plot=True)

    go_enr.results.to_csv(basefile + '.go_process.csv')
    shutil.rmtree(osp.join(basefile))

    go_enr = gp.enrichr(gene_list,
                          'GO_Cellular_Component_2018',
                          'human',
                          'go_output2',
                          osp.join(basefile, 'go_output'),
                          background,
                          no_plot=True)

    go_enr.results.to_csv(basefile + '.go_component.csv')
    shutil.rmtree(osp.join(basefile))

    go_enr = gp.enrichr(gene_list,
                          'GO_Molecular_Function_2018',
                          'human',
                          'go_output3',
                          osp.join(basefile, 'go_output'),
                          background,
                          no_plot=True)

    go_enr.results.to_csv(basefile + '.go_function.csv')
    shutil.rmtree(osp.join(basefile))

    wiki_enr = gp.enrichr(gene_list,
                        'WikiPathways_2019_Human',
                        'human',
                        'wiki_output',
                        osp.join(basefile, 'wiki_output'),
                        background,
                        no_plot=True)

    wiki_enr.results.to_csv(basefile + '.wikipathways.csv')
    shutil.rmtree(osp.join(basefile))

    reactome_enr = gp.enrichr(gene_list,
                        'Reactome_2016',
                        'human',
                        'wiki_output',
                        osp.join(basefile, 'reactome_output'),
                        background,
                        no_plot=True)

    reactome_enr.results.to_csv(basefile + '.reactome.csv')
    shutil.rmtree(osp.join(basefile))

    hpm_enr = gp.enrichr(gene_list,
                        'Tissue_Protein_Expression_from_Human_Proteome_Map',
                        'human',
                        'wiki_output',
                        osp.join(basefile, 'hpm_map_output'),
                        background,
                        no_plot=True)

    hpm_enr.results.to_csv(basefile + '.human_proteome_map.csv')
    shutil.rmtree(osp.join(basefile))

    proteomicsdb_enr = gp.enrichr(gene_list,
                        'Tissue_Protein_Expression_from_ProteomicsDB',
                        'human',
                        'wiki_output',
                        osp.join(basefile, 'proteomics_db_output'),
                        background,
                        no_plot=True)

    proteomicsdb_enr.results.to_csv(basefile + '.proteomics_db.csv')
    shutil.rmtree(osp.join(basefile))

    jensen_tissues_enr = gp.enrichr(gene_list,
                              'Jensen_TISSUES',
                              'human',
                              'wiki_output',
                              osp.join(basefile, 'jensen_tissues_output'),
                              background,
                              no_plot=True)

    jensen_tissues_enr.results.to_csv(basefile + '.jensen_tissues.csv')
    shutil.rmtree(osp.join(basefile))

    covid19_enr = gp.enrichr(gene_list,
                              'COVID-19_Related_Gene_Sets',
                              'human',
                              'wiki_output',
                              osp.join(basefile, 'covid19_output'),
                              background,
                              no_plot=True)

    covid19_enr.results.to_csv(basefile + '.covid19.csv')
    shutil.rmtree(osp.join(basefile))


class ErcWorkspace:

    def __init__(self, directory: str, timetree: str, segmented: bool = False, segment_size: int = 10,
                 internal_requirement: float = -1, include_terminal: bool = False, recalculate: bool = False,
                 sliding_window: bool = False, skip_align: bool = False, skip_trim: bool = False,
                 taxon_set: List[str] = None, time_corrected: bool = False, id2name: Dict[str, str] = dict()):
        self.directory = directory
        self.timetree = timetree
        self.aligns = []
        self.paired_aligns = []
        self.concatenated = dict()
        self.segmented = segmented
        self.segment_size = segment_size
        self.internal_requirement = internal_requirement
        self.include_terminal = include_terminal
        self.recalculate = recalculate
        self.sliding_window = sliding_window
        self.skip_align = skip_align
        self.skip_trim = skip_trim
        self.taxon_set = taxon_set
        self.time_corrected = time_corrected
        self.id2name = id2name

        safe_mkdir(directory)
        safe_mkdir(osp.join(directory, 'cleaned'))
        safe_mkdir(osp.join(directory, 'aligns'))
        safe_mkdir(osp.join(directory, 'trim'))
        safe_mkdir(osp.join(directory, 'pruned'))
        safe_mkdir(osp.join(directory, 'tree'))
        safe_mkdir(osp.join(directory, 'concat'))
        if segmented:
            safe_mkdir(osp.join(directory, 'segments'))

        os.chdir(directory)

    def add_alignment(self, alignment: str, alignment2: str = None):
        assert osp.exists(alignment)
        if not alignment2:
            self.aligns.append(alignment)
        else:
            assert osp.exists(alignment2)
            self.paired_aligns.append((alignment, alignment2))

    def add_concatenated_alignment(self, alignments: List[str], name: str,
                                   alignments2: List[str] = None, name2: str = None):
        if '.' not in name:
            name += ".fa"
        if name2 and '.' not in name2:
            name2 += ".fa"

        name = osp.basename(name)
        self.concatenated[name] = (alignments, osp.join(self.directory, 'concat', name + '.part'))
        if not alignments2:
            self.aligns.append(name)
        else:
            name2 = osp.basename(name2)
            self.concatenated[name2] = (alignments2, osp.join(self.directory, 'concat', name2 + '.part'))
            self.paired_aligns.append((name, name2))

    async def segment(self, alignment: str):
        input_path = osp.join(self.directory, 'trim', alignment)
        output_prefix = osp.join(self.directory, 'segments', alignment)
        records = read_records(input_path)
        length = max([len(x.sequence) for x in records])
        chunk_count = int(math.ceil(length / self.segment_size)) if not self.sliding_window else (length + 1 - self.segment_size)
        for c in range(chunk_count):
            chunk_records = []
            is_last = c == chunk_count - 1
            start_index = (c * self.segment_size) if not self.sliding_window else c
            segment_len = 0
            for rec in records:
                if start_index >= len(rec.sequence):
                    continue
                if is_last:
                    seq = rec.sequence[start_index:]
                else:
                    end = ((c+1)*self.segment_size) if not self.sliding_window else (start_index+self.segment_size)
                    seq = rec.sequence[start_index:end]

                non_gap_seq = "".join(char for char in seq if char not in GAP_CHARS)
                if len(non_gap_seq) <= 1:
                    continue

                segment_len = len(seq)

                chunk_records.append(Record(rec.title, seq))
            if len(chunk_records) <= 1:
                print(f"WARNING! Chunk {start_index}-{start_index+self.segment_size} from {alignment} is uninformative and being skipped!")
                continue
            write_records(output_prefix + f'.chunk_{start_index}_{min(self.segment_size, segment_len) + start_index}',
                          chunk_records)

    async def run(self):
        align_names = [osp.basename(f) for f in self.aligns]
        concat_align_names = []

        for pair in self.paired_aligns:
            if osp.basename(pair[0]) not in align_names:
                align_names.append(osp.basename(pair[0]))
            if osp.basename(pair[1]) not in align_names:
                align_names.append(osp.basename(pair[1]))

        if not osp.exists(osp.join(self.directory, 'aligns.tar.bz2')):
            for aligns in chunks(self.aligns, 4):
                await asyncio.gather(*[clean_fasta(self.timetree, f, osp.join(self.directory, 'cleaned', osp.basename(f))) for f in aligns if f not in self.concatenated])
            for aligns in chunks(self.paired_aligns, 4):
                await asyncio.gather(*[clean_fasta(self.timetree, f[0], osp.join(self.directory, 'cleaned', osp.basename(f[0]))) for f in aligns if f[0] not in self.concatenated])
                await asyncio.gather(*[clean_fasta(self.timetree, f[1], osp.join(self.directory, 'cleaned', osp.basename(f[1]))) for f in aligns if f[1] not in self.concatenated])

            # Add missing concatenated
            for key, val in self.concatenated.items():
                aligns = [a for a in val[0] if osp.basename(a) not in align_names and osp.basename(a) not in concat_align_names]
                concat_align_names += [osp.basename(a) for a in aligns]
                for als in chunks(aligns, 4):
                    await asyncio.gather(*[clean_fasta(self.timetree, f, osp.join(self.directory, 'cleaned', osp.basename(f))) for f in als])

        if not osp.exists(osp.join(self.directory, 'aligns.tar.bz2')):
            if self.skip_align:
                for f in (align_names + concat_align_names):
                    if f not in self.concatenated:
                        shutil.copy(osp.join(self.directory, 'cleaned', f), osp.join(self.directory, 'aligns', f))
            else:
                for aligns in chunks(align_names + concat_align_names, 4):
                    await asyncio.gather(*[align(osp.join(self.directory, 'cleaned', f),
                                                 osp.join(self.directory, 'aligns', f)) for f in aligns if f not in self.concatenated])
            shutil.rmtree(osp.join(self.directory, 'cleaned'))

        to_remove = set()  # Remove seqs only composed of gaps
        if not osp.exists(osp.join(self.directory, 'trim.tar.bz2')):
            for aligns in chunks(align_names + concat_align_names, 4):
                if self.skip_trim:
                    for f in aligns:  # Fake trimming
                        if f not in self.concatenated:
                            shutil.copy(osp.join(self.directory, 'aligns', f), osp.join(self.directory, 'trim', f))
                else:
                    await asyncio.gather(*[trim(osp.join(self.directory, 'aligns', f),
                                                osp.join(self.directory, 'trim', f)) for f in aligns if f not in self.concatenated])
                for f in aligns:
                    if not osp.exists(osp.join(self.directory, 'trim', f)):
                        if f not in self.concatenated:
                            to_remove.add(f)
            for f in self.concatenated.keys():
                if osp.exists(osp.join(self.directory, 'trim', f)):
                    continue
                concat_entry = self.concatenated[f]
                print("Concatenation: ", f, concat_entry)
                await concat([osp.join(self.directory, 'trim', osp.basename(e)) for e in concat_entry[0]],
                             concat_entry[1], osp.join(self.directory, 'trim', f))

            await archive_directory(osp.join(self.directory, 'aligns'), "aligns.tar.bz2")

        align_names = [a for a in align_names if a not in to_remove]

        if self.segmented:  # Note: Concatentation not considered in segments because its pointless
            if len(self.concatenated) > 0:
                raise AssertionError("Concatenation not supported with segERCs!")
            if not osp.exists(osp.join(self.directory, 'pruned.tar.bz2')):
                for aligns in chunks(align_names, 4):
                    await asyncio.gather(*[self.segment(a) for a in aligns])
                    for align_file in aligns:
                        for segments in chunks([f for f in os.listdir(osp.join(self.directory, 'segments')) if align_file in f], 4):
                            await asyncio.gather(*[prune(osp.join(self.directory, 'segments', f), self.timetree,
                                                         osp.join(self.directory, 'pruned', f)) for f in segments])

            for align_file in align_names:
                for segments in chunks([f for f in os.listdir(osp.join(self.directory, 'segments')) if align_file in f], 8):
                    await asyncio.gather(*[protein_tree(osp.join(self.directory, 'segments', f),
                                                        osp.join(self.directory, 'pruned', f),
                                                        osp.join(self.directory, 'tree', f + '.pred')) for f in segments if not osp.exists(osp.join(self.directory, 'tree', f + '.pred'))])
        else:
            if not osp.exists(osp.join(self.directory, 'pruned.tar.bz2')):
                for aligns in chunks(align_names, 4):
                    await asyncio.gather(*[prune(osp.join(self.directory, 'trim', f), self.timetree,
                                                 osp.join(self.directory, 'pruned', f)) for f in aligns])

            for aligns in chunks([n for n in align_names if n not in self.concatenated], 8):
                await asyncio.gather(*[protein_tree(osp.join(self.directory, 'trim', f),
                                                    osp.join(self.directory, 'pruned', f),
                                                    osp.join(self.directory, 'tree', f + '.pred')) for f in aligns if not osp.exists(osp.join(self.directory, 'tree', f + '.pred'))])
            for aligns in chunks(list(self.concatenated.keys()), 8):
                await asyncio.gather(*[protein_tree(osp.join(self.directory, 'trim', f),
                                                    osp.join(self.directory, 'pruned', f),
                                                    osp.join(self.directory, 'tree', f + '.pred'),
                                                    self.concatenated[f][1]) for f in aligns if not osp.exists(osp.join(self.directory, 'tree', f + '.pred'))])

        await archive_directory(osp.join(self.directory, 'trim'), 'trim.tar.bz2')
        await archive_directory(osp.join(self.directory, 'pruned'), 'pruned.tar.bz2')
        if self.segmented:
            await archive_directory(osp.join(self.directory, 'segments'), 'segments.tar.bz2')

        past_trees = set()

        for source in datasources:
            for entry in read_datasource(source):
                past_trees.add(entry.path1)
                past_trees.add(entry.path2)

        for past_tree in past_trees:
            shutil.copy(past_tree, osp.join(self.directory, 'tree', osp.basename(past_tree)))

        past_trees = {osp.basename(f) for f in past_trees}
        all_trees = {osp.basename(f) + ".pred" for f in self.aligns} | past_trees

        completed = list()
        if not osp.exists(osp.join(self.directory, 'ercs_raw.csv')):
            with open(osp.join(self.directory, 'ercs_raw.csv'), 'w') as f:
                f.write('file1,file2,taxon,file1_terminal,file2_terminal,time\n')
        if not osp.exists(osp.join(self.directory, 'ercs.csv')):
            with open(osp.join(self.directory, 'ercs.csv'), 'w') as f:
                f.write('file1,file2,erc,p\n')
        else:
            with open(osp.join(self.directory, "ercs.csv"), 'r') as f:
                first = True
                for l in f:
                    if first:
                        first = False
                        continue
                    try:
                        split = l.split(",")[:2]
                        completed.append((split[0], split[1]))
                    except:
                        pass

        if self.segmented:
            tree_combos = []
            for combo in itertools.combinations(all_trees, 2):
                chunked_trees = list(glob(osp.join(self.directory, 'tree', combo[0].replace(".pred", ".chunk*.pred")))) \
                                + list(glob(osp.join(self.directory, 'tree', combo[1].replace(".pred", ".chunk*.pred"))))
                for chunk1, chunk2 in itertools.combinations(chunked_trees, 2):
                    if (combo[0] not in past_trees or combo[1] not in past_trees or self.recalculate) and (combo[0] != combo[1]) \
                            and (combo not in completed) and ((combo[1], combo[0]) not in completed):
                        tree_combos.append((osp.basename(chunk1), osp.basename(chunk2)))
        else:
            tree_combos = [combo for combo in itertools.combinations(all_trees, 2)
                           if (combo[0] not in past_trees or combo[1] not in past_trees or self.recalculate) and (combo[0] != combo[1])
                           and (combo not in completed) and ((combo[1], combo[0]) not in completed)]

        tree_combos = [combo for combo in tree_combos if combo[0].split(".chunk")[0] != combo[1].split(".chunk")[0]]
        tree_combos = [combo for combo in tree_combos if osp.basename(combo[0]) not in to_remove and osp.basename(combo[1]) not in to_remove]

        for pair in self.paired_aligns:
            tree_pair = (osp.basename(pair[0]) + '.pred', osp.basename(pair[1]) + '.pred')
            if self.segmented:
                for chunk1 in glob(osp.join(self.directory, 'tree', tree_pair[0].replace(".pred", ".chunk*.pred"))):
                    for chunk2 in glob(osp.join(self.directory, 'tree', tree_pair[1].replace(".pred", ".chunk*.pred"))):
                        chunked_pair = (osp.basename(chunk1), osp.basename(chunk2))
                        tree_combos.append(chunked_pair)
            else:
                if len([combo for combo in tree_combos if combo[0] in tree_pair and combo[1] in tree_pair]) == 0:
                    tree_combos.append(tree_pair)

        del past_trees
        del all_trees
        del completed

        for filepairs in chunks(tree_combos, 4):
            results = await asyncio.gather(*[generate_erc(self.timetree, osp.join(self.directory, 'tree', f1),
                                                          osp.join(self.directory, 'tree', f2),
                                                          self.internal_requirement,
                                                          self.include_terminal,
                                                          self.taxon_set,
                                                          self.time_corrected) for (f1, f2) in filepairs])

            with open(osp.join(self.directory, 'ercs.csv'), 'a') as f:
                for res in results:
                    f.write(f"{osp.basename(res.file1)},{osp.basename(res.file2)},{res.rho},{res.p}\n")
            with open(osp.join(self.directory, 'ercs_raw.csv'), 'a') as f:
                for res in results:
                    if res.raw_rates is not None:
                        for taxon, raw in res.raw_rates.items():
                            f.write(f"{osp.basename(res.file1)},{osp.basename(res.file2)},{taxon},{raw[0]},{raw[1]},{raw[2]}\n")

        if self.segmented:
            segmented_erc_results = []
            safe_mkdir(osp.join(self.directory, 'segmented_ercs'))
            with open(osp.join(self.directory, 'ercs.csv'), 'r') as f:
                first = True
                for l in f:
                    if first:
                        first = False
                        continue
                    if len(l.strip()) == 0:
                        continue
                    split = l.split(",")
                    if split[0] < split[1]:
                        file1 = split[0]
                        file2 = split[1]
                    else:
                        file1 = split[1]
                        file2 = split[0]
                    res = ErcResult(file1, file2, np.float(split[2]), np.float(split[3]), None)
                    file1_chunked = '.chunk_' in res.file1
                    file2_chunked = '.chunk_' in res.file2

                    if not file1_chunked and not file2_chunked:
                        continue

                    if file1_chunked or file2_chunked:
                        file1_base = res.file1
                        file2_base = res.file2
                        file1_range = "N/A"
                        file2_range = "N/A"

                        if file1_chunked:
                            split = file1_base.split(".chunk_")
                            file1_base = split[0]
                            file1_range = split[1].split(".")[0]

                        if file2_chunked:
                            split = file2_base.split(".chunk_")
                            file2_base = split[0]
                            file2_range = split[1].split(".")[0]

                        segmented_erc_results.append(SegmentedERC(res, file1_chunked, file2_chunked,
                                                                  file1_base, file2_base, file1_range, file2_range))

            grouped = dict()
            for res in segmented_erc_results:
                key = (res.file1_base, res.file2_base)
                if key not in grouped:
                    grouped[key] = []
                grouped[key].append(res)

            for (input1, input2) in grouped.keys():
                res_name = f"{osp.basename(input1).split('.')[0]}_vs_{osp.basename(input2).split('.')[0]}"
                values = grouped[(input1, input2)]
                values: List[SegmentedERC] = list(sorted(values,
                                     key=lambda res: (int(res.file1_range.split("_")[0]) if res.is_file1_segmented else 0,
                                                      int(res.file2_range.split("_")[0]) if res.is_file2_segmented else 0)))

                with open(osp.join(self.directory, 'segmented_ercs', res_name + ".flat.tsv"), 'w') as f:
                    f.write("file1\tfile2\tfile1_range\tfile2_range\trho\tp\n")
                    for value in values:
                        f.write(f"{value.file1_base}\t{value.file2_base}\t{value.file1_range}\t{value.file2_range}\t{value.result.rho}\t{value.result.p}\n")

                connections = nx.Graph()
                file1_base = None
                file2_base = None
                for value in values:
                    if self.sliding_window and value.is_file1_segmented:
                        file1_split = value.file1_range.split("_")
                        for file1_i in range(int(file1_split[0]), int(file1_split[1])):
                            file1_i = "file1_" + str(file1_i)
                            if not connections.has_node(file1_i):
                                connections.add_node(file1_i)
                    else:
                        file1_range = "file1_" + value.file1_range
                        if not connections.has_node(file1_range):
                            connections.add_node(file1_range)

                    if self.sliding_window and value.is_file2_segmented:
                        file2_split = value.file2_range.split("_")
                        for file2_i in range(int(file2_split[0]), int(file2_split[1])):
                            file2_i = "file2_" + str(file2_i)
                            if not connections.has_node(file2_i):
                                connections.add_node(file2_i)
                    else:
                        file2_range = "file2_" + value.file2_range
                        if not connections.has_node(file2_range):
                            connections.add_node(file2_range)

                    if self.sliding_window:
                        file1_positions = ["file1_" + value.file1_range]
                        file2_positions = ["file2_" + value.file2_range]
                        if value.is_file1_segmented:
                            file1_positions = []
                            file1_split = value.file1_range.split("_")
                            for file1_i in range(int(file1_split[0]), int(file1_split[1])):
                                file1_i = "file1_" + str(file1_i)
                                file1_positions.append(file1_i)
                        if value.is_file2_segmented:
                            file2_positions = []
                            file2_split = value.file2_range.split("_")
                            for file2_i in range(int(file2_split[0]), int(file2_split[1])):
                                file2_i = "file2_" + str(file2_i)
                                file2_positions.append(file2_i)

                        for file1_i in file1_positions:
                            for file2_i in file2_positions:
                                if not connections.has_edge(file1_i, file2_i):
                                    connections.add_edge(file1_i, file2_i, rho=0, p=0, count=0)
                                connections[file1_i][file2_i]['count'] += 1
                                count = connections[file1_i][file2_i]['count']
                                connections[file1_i][file2_i]['rho'] += ((value.result.rho - connections[file1_i][file2_i]['rho'])/count)
                                connections[file1_i][file2_i]['p'] += ((value.result.p - connections[file1_i][file2_i]['p'])/count)
                    else:
                        file1_range = "file1_" + value.file1_range
                        file2_range = "file2_" + value.file2_range
                        connections.add_edge(file1_range, file2_range, rho=value.result.rho, p=value.result.p, count=1)

                    file1_base = value.file1_base
                    file2_base = value.file2_base

                file1_cols = list(sorted([n for n in connections.nodes if "file1_" in n], key=lambda n: int(n.split("_")[1])))
                file2_cols = list(sorted([n for n in connections.nodes if "file2_" in n], key=lambda n: int(n.split("_")[1])))
                with open(osp.join(self.directory, 'segmented_ercs', res_name + ".mat.p.tsv"), 'w') as f:
                    f.write(f"\"{file1_base}:\n{file2_base} V\"\t")
                    f.write("\t".join(file1_cols))
                    f.write("\n")
                    for file2_col in file2_cols:
                        f.write(f"{file2_col}")
                        for file1_col in file1_cols:
                            f.write(f"\t{connections[file1_col][file2_col]['p']}")
                        f.write("\n")
                with open(osp.join(self.directory, 'segmented_ercs', res_name + ".mat.rho.tsv"), 'w') as f:
                    f.write(f"\"{file1_base}:\n{file2_base} V\"\t")
                    f.write("\t".join(file1_cols))
                    f.write("\n")
                    for file2_col in file2_cols:
                        f.write(f"{file2_col}")
                        for file1_col in file1_cols:
                            f.write(f"\t{connections[file1_col][file2_col]['rho']}")
                        f.write("\n")
                with open(osp.join(self.directory, 'segmented_ercs', res_name + ".mat.padjsig.rho.tsv"), 'w') as f:
                    f.write(f"\"{file1_base}:\n{file2_base} V\"\t")
                    f.write("\t".join(file1_cols))
                    f.write("\n")
                    for file2_col in file2_cols:
                        f.write(f"{file2_col}")
                        for file1_col in file1_cols:
                            if np.float(connections[file1_col][file2_col]['p']) * connections.number_of_edges() < 0.05:
                                f.write(f"\t{connections[file1_col][file2_col]['rho']}")
                            else:
                                f.write("\t0.0")
                        f.write("\n")


    def generate_network(self) -> nx.Graph:
        """Assumes that we are still in the workspace directory"""
        graph = nx.Graph()
        for entry in read_datasource(ErcDataSource("tree/", "ercs.csv")):
            id1 = osp.basename(entry.path1).split('.')[0]
            id2 = osp.basename(entry.path2).split('.')[0]

            graph.add_edge(id1, id2, rho=entry.rho, p=entry.p, weight=entry.rho)
        return graph


    def export_results(self, basename: str):
        book = xlsxwriter.Workbook(basename)
        rho_sheet = book.add_worksheet("Rho Values")
        p_sheet = book.add_worksheet("P Values")

        bold = make_bold_formatting(book)
        p_format = make_p_formatting(book)
        rho_format = make_rho_formatting(book)

        network = self.generate_network()

        sorted_ids = list(sorted(network.nodes, key = lambda n: self.id2name.get(n, n)))

        rho_sheet.write_row(0, 1, [self.id2name.get(n, n) for n in sorted_ids], cell_format=bold)
        p_sheet.write_row(0, 1, [self.id2name.get(n, n) for n in sorted_ids], cell_format=bold)
        rho_sheet.write_column(1, 0, [self.id2name.get(n, n) for n in sorted_ids], cell_format=bold)
        p_sheet.write_column(1, 0, [self.id2name.get(n, n) for n in sorted_ids], cell_format=bold)

        for i, id1 in enumerate(sorted_ids):
            for j, id2 in enumerate(sorted_ids):
                if id1 == id2:
                    rho_sheet.write(i + 1, j + 1, 1.0)
                    p_sheet.write(i + 1, j + 1, 1.0)
                else:
                    rho_sheet.write(i + 1, j + 1, network[id1][id2]['rho'], rho_format)
                    p_sheet.write(i + 1, j + 1, network[id1][id2]['p'], p_format)

        book.close()
