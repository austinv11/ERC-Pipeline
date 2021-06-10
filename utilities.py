import asyncio
import aiohttp
import itertools
import math
import os
import os.path as osp
import sys
import traceback
from typing import Tuple, Dict, List, Callable, Any, Union

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


def decorate(fun: Callable[..., Any]) -> type:
    class decorator:

        def __init__(self):
            pass

        def __call__(self, f):
            def wrapped(*args, **kwargs):
                try:
                    o = fun(*args, **kwargs)
                    if o is not None:
                        return o
                except:
                    pass
                return f(*args, **kwargs)

            return wrapped

    return decorator


RECURSION_LIMIT = 5


class Limiter:

    def __init__(self, delay: int = 1):
        self._time = -1
        self._limiter_lock = asyncio.Lock()
        self.request_delay = delay
        self.changed = False

    async def __aenter__(self):
        await self._limiter_lock.acquire()

        to_wait = self._time_to_wait()

        await asyncio.sleep(to_wait)

    async def __aexit__(self, *args, **kwargs):
        self._update_timer()

        self._limiter_lock.release()

    def _time_to_wait(self) -> int:
        request_time = self._time

        if request_time == -1:
            return 0

        now = asyncio.get_event_loop().time()

        to_wait = request_time + self.request_delay - now
        to_wait = max(0, to_wait)

        return to_wait

    def _update_timer(self):
        now = asyncio.get_event_loop().time()
        self._time = now


_new_session_lock = asyncio.Lock()
_limiters = dict()

_sessions = dict()


async def remote_call(url: str, is_json: bool, *, recursion_count: int = 0, **kwargs: Any) -> Union[Dict, str]:
    orig_url = url
    if recursion_count >= RECURSION_LIMIT:
        raise Exception(f"Too much recursion! (url={url})")

    if kwargs is not None and len(kwargs) > 0:
        first = True
        for k, v in kwargs.items():
            if v is None:
                continue

            if first:
                url = url + '?'
                first = False
            else:
                url = url + '&'

            if isinstance(v, set) or isinstance(v, list):
                v = ",".join(v)

            url = url + f'{k}={v}'

    base_url = url.split('/')[0]

    await _new_session_lock.acquire()

    if base_url not in _sessions:
        _sessions[base_url] = aiohttp.ClientSession()
        _limiters[base_url] = Limiter()

    session = _sessions[base_url]
    limiter = _limiters[base_url]

    _new_session_lock.release()

    async with limiter:
        try:
            async with session.get(url) as resp:
                resp.raise_for_status()
                if not limiter.changed and resp.headers.get('x-rate-limit-limit', None) is not None:  # Prepare ratelimiting if applicable
                    req_limit = float(resp.headers.get('x-rate-limit-limit'))  # We are assuming a 60 second window
                    limit_time = float(resp.headers.get('x-rate-limit-reset', "1"))  # We are assuming that the first request will give us the window time
                    limiter.request_delay = limit_time / req_limit
                    limiter.changed = True
                if resp.status != 503:
                    if resp.status == 429:  # Rate limit
                        if int(resp.headers.get('x-rate-limit-remaining', 0)) == 0:
                            await asyncio.sleep(float(resp.headers.get('x-rate-limit-reset', "1")))
                    elif is_json:
                        return await resp.json(content_type=None)
                    else:
                        return await resp.text()
        except Exception as e:
            print(f"Exception for url {url} caught!", flush=True)
            traceback.print_exc()
            await asyncio.sleep(1)  # Sleep a little

    return await remote_call(orig_url, is_json, recursion_count=recursion_count + 1, **kwargs)


async def close_web_sessions():
    for session in _sessions.values():
        await session.close()
    _sessions.clear()
