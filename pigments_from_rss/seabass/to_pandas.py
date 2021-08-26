import pandas as pd

from .reader import readSB


def seabass_to_pandas(path):
    """SeaBASS to Pandas DataFrame converter
    
    Parameters
    ----------
    path : str
        path to an FCHECKed SeaBASS file
    
    Returns
    -------
    pandas.DataFrame
    
    """
    sb = readSB(path)
    dataframe = pd.DataFrame.from_dict(sb.data)
    return dataframe
