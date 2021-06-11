import os.path as osp
import sys

from statsmodels.stats.multitest import multipletests

from pipeline import register_previous_run, generate_network, _30mya_cutoff, _20mya_cutoff, find_tree, get_rates
from utilities import safe_mkdir, rho_sorted_neighbors, safe_phylo_read, mammal_taxa_info


def main():
    args = sys.argv[1:]

    assert len(args) >= 5, "Arguments: python3 get_rate_data.py full/20my/30my pipeline_directory output protein_id [protein_id 2]"

    tree = safe_phylo_read(osp.join("data", "finished_mam_timetree.nwk"))
    taxon_set = args[0].lower()

    assert taxon_set in ("full", "20my", "30my"), "You must specify the taxon set to get rates for."

    taxa_list = None
    if taxon_set == '20my':
        taxa_list = _20mya_cutoff
    elif taxon_set == '30my':
        taxa_list = _30mya_cutoff

    pipeline_dir = args[1]
    assert osp.exists(pipeline_dir), f"Pipeline directory {pipeline_dir} doesn't exist!"
    register_previous_run(pipeline_dir)

    output = args[2]
    assert output.endswith(".csv"), f"The output ({output}) must be a .csv file!"

    protein = args[3].split(".")[0].strip()

    if len(args) > 5:
        protein2 = args[4].split(".")[0].strip()
    else:
        protein2 = None

    protein = find_tree(protein)
    if protein2:
        protein2 = find_tree(protein2)

    name2taxa = mammal_taxa_info(name_as_key=True)

    if protein2:
        print(f"Getting rate data for {protein} with {protein2}...")
        taxa, rates = get_rates(tree, True, taxa_list, protein, protein2)
        with open(output, 'w') as f:
            f.write("protein1,protein2,taxon,taxon order,rate1,rate2,time\n")
            for (taxon, rate1, rate2, time) in zip(taxa, rates[1], rates[2], rates[0]):
                order = name2taxa[taxon.strip()].order
                f.write(f"{protein},{protein2},{taxon},{order},{rate1},{rate2},{time}\n")
    else:
        print(f"Getting rate data for {protein}...")
        taxa, rates = get_rates(tree, True, taxa_list, protein)
        with open(output, 'w') as f:
            f.write("protein,taxon,taxon order,time,rate\n")
            for (taxon, rate1, time) in zip(taxa, rates[1], rates[0]):
                order = name2taxa[taxon.strip()].order
                f.write(f"{protein},{taxon},{order},{time},{rate1}\n")


if __name__ == "__main__":
    main()
