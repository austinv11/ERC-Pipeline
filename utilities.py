import asyncio
import itertools
import math
import os
import os.path as osp
import sys
from typing import Tuple, Dict, List

import xlsxwriter
from ete3 import PhyloTree
from ete3.parser.newick import NewickError
import networkx as nx


def safe_phylo_read(filename) -> PhyloTree:
    try:
        return PhyloTree(filename, format=3)
    except:
        try:
            return PhyloTree(filename)
        except:
            try:
                return PhyloTree(filename, format=1)
            except:
                try:
                    return PhyloTree(filename, format=5)
                except NewickError as e:
                    print(f"Are you sure tree {filename} exists?", file=sys.stderr, flush=True)
                    raise e


def safe_delete(path: str):
    if osp.exists(path):
        try:
            os.remove(path)
        except: pass


def safe_mkdir(dir: str):
    if not osp.exists(dir):
        try:
            os.makedirs(dir)
        except: pass


def wait(coro):  # Wait on coroutines more simply
    return asyncio.get_event_loop().run_until_complete(coro)


# noinspection PyPep8Naming
def translate(G: nx.Graph, l2n: Dict[str, str]) -> nx.Graph:
    G_c = nx.Graph() if not G.is_directed() else nx.DiGraph()

    for u, v in G.edges:
        G_c.add_edge(l2n.get(u, u), l2n.get(v, v), **G.get_edge_data(u, v), default=dict())
    return G_c


def override_sys_out(tag: str = None):
    import sys
    from datetime import datetime as dt

    def new_write(iostream):

        orig = iostream.write

        class TimeStamper:

            nl = True

            def write(self, s):
                """Write function overloaded."""
                if s == '\n':
                    orig(s)
                    self.nl = True
                elif self.nl:
                    orig('[%s]%s: %s' % (str(dt.now()), '' if not tag else f"[{tag}]", s))
                    self.nl = False
                else:
                    orig(s)

        stamper = TimeStamper()

        iostream.write = stamper.write

    new_write(sys.stdout)
    new_write(sys.stderr)


async def async_call(cmd, cwd=None, **kwargs):
    print("> " + cmd, flush=True)
    if "shell" in kwargs.keys():
        del kwargs["shell"]
    kwargs['cwd'] = cwd if cwd else os.getcwd()

    process = await asyncio.create_subprocess_shell(cmd, **kwargs)
    await process.wait()


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


def try_hook_uvloop():  # Can improve speed on unix systems
    try:
        import uvloop
        uvloop.install()
    except: pass


def _self_path():
    path = osp.dirname(__file__)
    if not path:
        path = '.'
    return path


def _add_edge(net: nx.Graph, a, b):
    if not net.has_edge(a, b):
        net.add_edge(a, b)


def rho_sorted_neighbors(net: nx.Graph, node: str) -> List[str]:
    return list(sorted(net.neighbors(node), key=lambda n: net[n][node]['rho'], reverse=True))


def p_sorted_neighbors(net: nx.Graph, node: str) -> List[str]:
    return list(sorted(net.neighbors(node), key=lambda n: net[n][node]['p']))


def rrn_nets(net: nx.Graph, *odbs: str) -> Tuple[nx.DiGraph, nx.DiGraph, nx.DiGraph]:
    step1 = nx.DiGraph()
    added_set = set()

    # Reciprank 20 first
    for prot in odbs:
        for i, neighbor in enumerate(rho_sorted_neighbors(net, prot)[:20]):
            prot_pos = rho_sorted_neighbors(net, neighbor).index(prot)
            if prot_pos < 20:
                added_set.add(neighbor)
                _add_edge(step1, prot, neighbor)
                _add_edge(step1, neighbor, prot)

    # Reciprank 20 of added
    step2 = step1.copy()
    for added in added_set:
        for i, neighbor in enumerate(rho_sorted_neighbors(net, added)[:20]):
            prot_pos = rho_sorted_neighbors(net, neighbor).index(added)
            if prot_pos < 20:
                _add_edge(step2, added, neighbor)
                _add_edge(step2, neighbor, added)

    # Fill connections <=20 between current node set
    step3 = step2.copy()
    for (n1, n2) in itertools.permutations(step2.nodes, 2):
        if n1 == n2 or step2.has_edge(n1, n2):
            continue
        prot_pos = rho_sorted_neighbors(net, n1).index(n2)
        if prot_pos < 20:
            _add_edge(step3, n1, n2)

    return step1, step2, step3


def make_rho_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        "num_format": "0.000"
    })


def make_p_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        "num_format": "0.00E+00"
    })


def make_bold_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        'bold': True
    })


def color_from_custom_map(pct: float, cmap: Dict[float, Tuple[int, int, int]]) -> str:
    # https://stackoverflow.com/a/7128796/5179044
    maxColor = 0
    lastMaxColor = 0
    for col in sorted(cmap.keys()):
        lastMaxColor = maxColor
        maxColor = col
        if pct < col:
            break

    _range = maxColor - lastMaxColor
    range_pct = (pct - lastMaxColor) / _range

    pct_lower = 1 - range_pct
    pct_upper = range_pct

    color = (
        math.floor(cmap[lastMaxColor][0] * pct_lower + cmap[maxColor][0] * pct_upper),
        math.floor(cmap[lastMaxColor][1] * pct_lower + cmap[maxColor][1] * pct_upper),
        math.floor(cmap[lastMaxColor][2] * pct_lower + cmap[maxColor][2] * pct_upper),
    )

    hex = "#{:02x}{:02x}{:02x}".format(color[0], color[1], color[2])

    return hex


def graphviz_network_plot(net: nx.Graph, output: str, highlight: Dict[str, str] = dict(), circo: bool = False,
                          highlight_by_font: bool = False):
    # output: a png file
    # Requires pygraphviz

    cmap = {0.0: (255, 255, 255),
            0.5: (204, 217, 255),
            1.0: (20, 139, 255)}

    max_connectivity = max([len(list(net.neighbors(n))) for n in net.nodes])

    A = nx.nx_agraph.to_agraph(net)
    A.graph_attr.update(strict=False, overlap=False, splines='true')
    A.node_attr['style'] = 'filled'
    if circo:
        A.layout(prog="circo")
    else:
        A.layout()
    for n in net.nodes:
        n2 = A.get_node(n)
        connectivity = len(list(net.neighbors(n)))
        n2.attr["fillcolor"] = color_from_custom_map(connectivity / max_connectivity, cmap) if highlight_by_font or (not highlight_by_font and n not in highlight) else highlight[n]
        if highlight_by_font and n in highlight:
            n2.attr['fontcolor'] = highlight[n]
    A.draw(output, args='-Gsize=10,10 -Gdpi=300')
    A.draw(output.replace(".png", ".svg"), format='svg')
