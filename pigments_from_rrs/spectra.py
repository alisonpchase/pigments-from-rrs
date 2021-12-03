from .to_pandas import seabass_to_pandas
import pandas as pd 

G1 = 0.0949  # g1 and g2 are values from Gordon et al., 1988
G2 = 0.0794
LNOT = 400  # reference lambda wavelength (nm)

def load_data(path):
    """Load data from file

    Parameters
    ----------
    path: a valid file path to the dataset

    Returns
    ------
    A pandas dataframe

    """
    return pd.read_csv(path)
    
