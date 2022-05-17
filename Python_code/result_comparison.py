# standard library dependencies
import os
import re 
import pip
import sys
import subprocess
from typing import Tuple, Mapping, List, Union

# external dependencies; already installed in Google colab environment
import pandas as pd 

# external; may require installing in Google colab environment
try:
    import matplotlib.pyplot as plt 
except ModuleNotFoundError:
    print("Installing matplotlib with pip...; if this doesn't work, run `!pip uninstall -y matplotlib` and then `!pip install matplotlib` directly in a Colab cell")
    pip_version = list(map(int, pip.__version__.split(".")))
    if pip_version[0] < 20:
        os.system('pip install matplotlib')
    else:
        pip.main(['install', '--user', 'matplotlib'])
    import matplotlib.pyplot as plt 

try:
    from supervenn import supervenn
except ModuleNotFoundError:
    print("Installing supervenn with pip...; if this doesn't work, run `!pip uninstall -y supervenn` and then `!pip install supervenn` directly in a Colab cell")
    pip_version = list(map(int, pip.__version__.split(".")))
    if pip_version[0] < 20:
        os.system('pip install supervenn')
    else:
        pip.main(['install', '--user', 'supervenn'])
    from supervenn import supervenn

# TODO: 1. add a filename filter for load_files_from_dir

def load_files_from_dir(dirpath:str, **read_csv_kwargs) -> Mapping[str, pd.DataFrame]:
    """
    Wrapper around pd.read_csv to load all files in `dirpath`
    using the pd.read_csv kwargs specified with `read_csv_kwargs`.
    ----------
    Arguments:
        dirpath:    string path to the directory whose files need to be loaded. 
                    should be the absolute path starting with `/content/drive...`
                    if running in Google colab.
        read_csv_kwargs:    kwargs for pd.read_csv. 
                            defaults to {"sep":'\t', "index_col":None} apart from
                            default settings (see References).
    ----------
    Returns:
        dataframes_dictionary: dictionary of (filename:pd.DataFrame) key:value pairs.
    ----------
    References:
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html?highlight=read_csv#pandas.read_csv
    """
    assert os.path.isdir(dirpath)
    files = os.listdir(dirpath)
    
    try:
        assert files != []
    except AssertionError as ae:
        message = "\n".join([
            "Could not find any files in: ",
            dirpath
        ])
        raise FileNotFoundError(message) from ae 

    read_csv_args = {
        "sep":'\t',
        "index_col":None
    }

    read_csv_args.update(read_csv_kwargs)

    dataframes_dictionary = {
        filename: pd.read_csv(
            os.path.join(dirpath, filename), 
            **read_csv_args
        )
        for filename in files
    }

    return dataframes_dictionary

def compare_with_max_adjp(  dataframes_dictionary:Mapping[str,  pd.DataFrame],
                            max_adjp:Union[List[float], float],
                            tool_names:List[str]=None,
                            statistic_column_names:List[str]=None):
    """
    Compares the set of differentially expressed genes at the confidence level 
    (FDR/adjusted p-value specified by `max_adjp`.
    ----------
    Arguments:
        dataframes_dictionary:  dictionary of (filename:pd.DataFrame) key:value 
                                pairs. 
        max_adjp:   float (or list of floats) representing the confidence level.
        tool_names: list of strings to use as tool labels. 
                    defaults to ["DESeq2", "edgeR"]
    ----------
    Returns:
        Nothing
    ----------
    References:
        https://pypi.org/project/supervenn/
    """
    if isinstance(max_adjp, float):
        max_adjp = [max_adjp]

    if tool_names is None:
        tool_names = ["DESeq2", "edgeR"]

    if statistic_column_names is None:
        statistic_column_names = ["padj", "edgeR.exactTest_FDR"]

    assert len(tool_names) == len(statistic_column_names)

    for alpha in max_adjp:
        all_filenames = list(dataframes_dictionary.keys())
        common_prefix = os.path.commonprefix(all_filenames)
        for filename, df in dataframes_dictionary.items():
            comparison = filename.replace(common_prefix,'')
            supervenn_args = []

            for result_columns, tool in zip(statistic_column_names, tool_names):
                degset = df.loc[df[result_columns] <= alpha, "Ensembl ID"]
                degset = set(degset.values.tolist())
                supervenn_args.append(
                    (f"{comparison} with {tool} @ FDR={alpha}", degset)
                )
            x = supervenn(
                [set_ for (name, set_) in supervenn_args ],
                set_annotations=[ name for (name, set_) in supervenn_args ]
            )
            plt.show()
