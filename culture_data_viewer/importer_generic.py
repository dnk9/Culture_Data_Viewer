import pandas as pd
from pathlib import Path
from typing import Dict


def extract(reactor_data_filename: Path) -> Dict[str, pd.DataFrame]:
    """
    Extract reactor data from a xlsx or csv file.

    In case of xlsx files each cdf will be named as the sheet it was taken from.
    As one csv file can contain only data for one culture, the resulting cdf_dict will have just one cdf, with key "1".

    Parameters
    ----------
    reactor_data_filename: Path
        The Path of the file to open

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with a cdf for each reactor.

    Raises
    ------
    RuntimeWarning
        If the file is not recognised.

    """
    extension = reactor_data_filename.suffix[1:]

    if extension == "csv":
        reactor_data_df = pd.read_csv(reactor_data_filename)
        cdf_dict = {"1": reactor_data_df}

    elif extension == "xlsx":
        cdf_dict = pd.read_excel(reactor_data_filename, sheet_name=None)

    else:
        raise TypeError("Input file not recognised. Should be 'xlsx' or 'csv'. ")

    return cdf_dict
