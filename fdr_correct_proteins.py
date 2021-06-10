import os.path as osp
import sys

from statsmodels.stats.multitest import multipletests

from pipeline import register_previous_run, generate_network
from utilities import safe_mkdir, rho_sorted_neighbors


def main():
    args = sys.argv[1:]

    assert len(args) > 2, "Arguments: python3 fdr_correct_proteins.py l2n.tsv output_directory pipeline1_directory pipeline2_directory..."

    assert osp.exists(args[0])
    id2name = dict()
    with open(args[0], 'r') as f:
        first = True
        for l in f:
            if first:
                first = False
                continue
            split = l.strip().split("\t")
            split = [s.strip().strip('"').strip("'") for s in split]
            if len(split) < 2:
                continue
            id2name[split[0].strip()] = split[1].strip()

    print(f"Correcting p-values and generating lists to {args[1]}...")

    safe_mkdir(args[1])

    for pipeline in args[2:]:
        register_previous_run(pipeline)

    net = generate_network()

    node2sorted = dict()
    for n in net.nodes:
        node2sorted[n] = rho_sorted_neighbors(net, n)

    for n in net.nodes:
        name = id2name.get(n, n)
        with open(osp.join(args[1], f"{name}_corrected_list.csv"), 'w') as f:
            f.write(f"ODB,Protein,Rank in {name}'s list,{name}'s rank in partner's list,rho,raw p,FDR p\n")
            fdr_corrected = multipletests([net[n][n2]['p'] for n2 in node2sorted[n]], method="fdr_bh")[1]
            for index1, n2 in enumerate(node2sorted[n]):
                index2 = node2sorted[n2].index(n)
                f.write(f"{n2},{id2name.get(n2, n2)},{index1 + 1},{index2 + 1},{net[n][n2]['rho']},{net[n][n2]['p']},{fdr_corrected[index1]}\n")

    print("Complete!")


if __name__ == "__main__":
    main()
