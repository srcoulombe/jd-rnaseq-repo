# standard library dependencies
import pip
pip_version = list(map(int, pip.__version__.split(".")))

import os 
from typing import List, Tuple, Set, Union

# external dependencies
import pandas as pd 

# external; may require installing in Google colab environment
try:
    from gprofiler import GProfiler
except ModuleNotFoundError:
    print("Installing gprofiler-official with pip...")
    if pip_version[0] < 20:
        os.system('pip install gprofiler-official')
    else:
        pip.main(['install', '--user', 'gprofiler-official'])
    from gprofiler import GProfiler


try:
    import networkx as nx
except ModuleNotFoundError:
    print("Installing networkx with pip...")
    if pip_version[0] < 20:
        os.system('pip install networkx==2.3')
    else:
        pip.main(['install', '--user', 'networkx==2.3'])
    import networkx as nx

try:
    import goenrich
except ModuleNotFoundError:
    print("Installing goenrich from srcoulombe's repo...")
    os.system('git clone https://github.com/srcoulombe/goenrich.git')
    os.system('pip install -e goenrich')
    import goenrich

try:
    import bokeh
except ModuleNotFoundError:
    print("Installing bokeh with pip...")
    if pip_version[0] < 20:
        os.system('pip install bokeh')
    else:
        pip.main(['install', '--user', 'bokeh'])
finally:
    from bokeh.io import show, output_file, output_notebook
    from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxSelectTool, BoxZoomTool, ResetTool, PanTool
    from bokeh.models.graphs import from_networkx, NodesAndLinkedEdges, EdgesAndLinkedNodes
    from bokeh.palettes import Spectral4

def run_gProfiler_gOST( gene_symbols_list:List[str], 
                        organism:str='hsapiens') -> pd.DataFrame:
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism=organism,
            query=gene_symbols_list)
    return df

def generate_go_tree(obo_file_path:str=None, gaf_file_path:str=None) -> nx.Graph:
    """
    """
    if obo_file_path is None: 
        obo_file_path = os.path.join(
            os.getcwd(), 
            'go-basic.obo'
        )
    onto = goenrich.obo.ontology(obo_file_path)

    if gaf_file_path is None:
        gaf_file_path = os.path.join(
            os.getcwd(), 
            'goa_human.gaf.gz'
        )
        
    annot = goenrich.read.goa(gaf_file_path)

    values = {k: set(v) for k,v in annot.groupby('go_id')['db_object_symbol']}

    # propagate the background through the ontology
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(onto, values, background_attribute)
    return onto

def plot_subgraph(  onto_graph:nx.Graph, node_subset:Union[List[str],Set[str]],
                    in_IPython:bool=False):
    """
    """
    G = onto_graph.subgraph(node_subset)
    # Show with Bokeh
    plot = Plot(plot_width=1600, plot_height=800,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))

    plot.title.text = "Graph Interaction Demonstration"

    node_hover_tool = HoverTool(tooltips=[("index", "@index"), ("club", "@club")])
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool(), TapTool(), BoxSelectTool(), PanTool())

    graph_renderer = from_networkx(G, nx.kamada_kawai_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
    graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
    graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])

    graph_renderer.edge_renderer.glyph = MultiLine(line_alpha=0.8, line_width=1)
    graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
    #graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

    graph_renderer.selection_policy = NodesAndLinkedEdges()
    #graph_renderer.inspection_policy = EdgesAndLinkedNodes()


    plot.renderers.append(graph_renderer)

    #output_file("interactive_graphs.html")
    if in_IPython: output_notebook()

    show(plot)