import asyncio
import os
import os.path as osp

import click

from utilities import override_sys_out, try_hook_uvloop, _self_path
from pipeline import ErcWorkspace, _20mya_cutoff, _30mya_cutoff


FASTA_ENDINGS = [".fa", ".faa", ".fs", ".fasta"]


@click.command()
@click.option("timetree", default=osp.join(_self_path(), "data", "finished_mam_timetree.nwk"),
              help="The time-scaled species phylogeny in newick format. Defaults to mammalian tree")
@click.option("alignments", default=".", help="The directory to search for alignments in.")
@click.option("wd", default=".", help="The working directory to write files to.")
@click.option("erc_type", type=click.Choice(['original', 'bt', "20my", "30my"], case_sensitive=False), default="original",
              help="The type of ERC to calculate. Note that the 20MY and 30MY ERC calculations expect the mammalian "
                   "data presented in Varela et al (2021).")
def main(timetree: str, alignments: str, wd: str, erc_type: str):
    try_hook_uvloop()
    override_sys_out("ERC")

    assert osp.exists(timetree)

    arg_mods = dict()
    if erc_type == "bt":
        arg_mods['time_corrected'] = True
    elif erc_type == "20my":
        arg_mods['taxon_set'] = _20mya_cutoff
    elif erc_type == "30my":
        arg_mods['taxon_set'] = _30mya_cutoff

    workspace = ErcWorkspace(wd, timetree, **arg_mods)

    for f in os.listdir(alignments):
        if osp.isdir(f):
            continue
        if f.split(".")[-1].lower() not in FASTA_ENDINGS:
            print(f"WARNING: Skipping alignment file {f}! It does not have a fasta file ending!")
            continue
        workspace.add_alignment(osp.join(alignments, f))

    asyncio.get_event_loop().run_until_complete(workspace.run())

    workspace.export_results("all2all_matrix")

    print("Completed!")


if __name__ == "__main__":
    main()
