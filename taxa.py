import os
import shutil

import requests

from utilities import decorate

from ete3 import NCBITaxa

DEFAULT_TAXADB = "./.etetoolkit/taxa.sqlite"
TAXA_DUMP = "./taxdump.tar.gz"

__is_setup__ = False


def setup():
    global __is_setup__
    if not __is_setup__:
        if not os.path.exists(DEFAULT_TAXADB):
            req = requests.get('http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz', stream=True)
            with open(TAXA_DUMP, 'wb') as f:
                shutil.copyfileobj(req.raw, f)
            NCBITaxa(dbfile=DEFAULT_TAXADB, taxdump_file=TAXA_DUMP)
        __is_setup__ = True


ncbi_taxa = decorate(setup)


@ncbi_taxa()
def get_taxa_info() -> NCBITaxa:
    return NCBITaxa(dbfile=DEFAULT_TAXADB, taxdump_file=None if os.path.exists(DEFAULT_TAXADB) else TAXA_DUMP)


# Warning: in order for this to work you must go into the venv folder and then enter the
# lib/python3.6/site-packages/ete3/ncbi_taxonomy/SQLite-Levenshtein sub directory and run "make"
@ncbi_taxa()
def scientific_name_to_txid(name):
    taxa = get_taxa_info()
    return taxa.get_fuzzy_name_translation(name)[0]


@ncbi_taxa()
def get_name_for_rank(txid, rank, default="no_name_found"):
    try:
        taxa = get_taxa_info()
        lineage2ranks = taxa.get_rank(taxa.get_lineage(txid))
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        new_txid = ranks2lineage.get(rank, default)
        translator = taxa.get_taxid_translator([new_txid])
        return translator[new_txid]
    except:
        return default


@ncbi_taxa()
def get_lineage_info(txid):
    return {
        "order": get_name_for_rank(txid, "order"),
        "suborder": get_name_for_rank(txid, "suborder"),
        "infraorder": get_name_for_rank(txid, "infraorder"),
        "parvorder": get_name_for_rank(txid, "parvorder"),
        "superfamily": get_name_for_rank(txid, "superfamily"),
        "family": get_name_for_rank(txid, "family"),
        "subfamily": get_name_for_rank(txid, "subfamily"),
        "genus": get_name_for_rank(txid, "genus"),
        "species": get_name_for_rank(txid, "species")
    }
