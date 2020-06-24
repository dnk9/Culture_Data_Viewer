import pandas as pd
import numpy as np
from openpyxl import load_workbook
from typing import Dict
import logging

logging.basicConfig(level=logging.DEBUG)


def open_offline_file(offln_file: str, sheet_name: str = None):
    """
    Opens the "offline_file" and outputs a dictionary of dataframes: one for each reactor/condition
    The "offline_file" is the "nice" excel file, in which data from OD and HPLC have already been merged.
    This file is very specific to my personal use case. It is constructed from the table coming from the hplc.
    Data coming from different reactors are all stored one below the other.
    """
    logging.debug(f'Running MERGE.open_offline_file with: \n offln_file: {offln_file} \n sheet_name: {sheet_name} ')

    # Checking if the xlsx file is fine: does it have the correct sheet?
    if sheet_name is None:
        # If no sheet name is specified, take just the first sheet in the work book
        offln_df = pd.read_excel(offln_file, sheet_name=0).dropna(subset=["Time (h)"])

    elif sheet_name in load_workbook(filename=offln_file).sheetnames:
        # This may be a bit redundant: opening the file twice
        offln_df = pd.read_excel(offln_file, sheet_name=sheet_name).dropna(subset=["Time (h)"])

    else:
        logging.warning(f'Opening of offline file failed: no sheet {sheet_name} in {offln_file}', exc_info=True)
        return None

    offln_df['InoculationTime []'] = pd.to_timedelta(offln_df["Time (h)"], unit="h")
    offln_df["Sample"] = True

    logging.debug(f'Reactors data found in the current offline file: {offln_df["Condition"].unique()}')
    offln_df_dict = dict(tuple(offln_df.groupby("Condition")))
    return offln_df_dict


def merge_cdf_to_offln_df(cdf: pd.DataFrame, offln_df: pd.DataFrame) -> pd.DataFrame:
    """
    merges a single cdf to a single offln_df. The join is a left join. Using merge_asof. Not adding new
    keys to the left dataframe.
    Implies that all the timings are correct and inside the left dataframe
    """
    logging.debug(f'Running function *MERGE.merge_cdf_to_offln_df*')

    merged = pd.merge_asof(cdf, offln_df, tolerance=pd.Timedelta('1min'), on='InoculationTime []')

    # Finding the duplicates
    duplicated = merged.duplicated(subset='Time (h)')

    # Removing the duplication, by setting back the offln_df data to NaN
    merged.loc[duplicated, offln_df.columns.drop('InoculationTime []')] = np.NaN

    # Setting the "Sample" status to False for each non-sample
    merged["Sample"].fillna(False, inplace=True)

    return merged


def merge_df_dicts(cdf_dict: Dict[str, pd.DataFrame], offln_df_dict: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Merge each cdf to the corresponding offln_df.

    Use "merge_cdf_to_offln_df" for each culture in the cdf_dict and offln_df_dict. The first six
    columns in the offln_df are ignored, as they represent information about the experimental conditions.
    """
    logging.debug(f'Running *MERGE.merge_df_dicts*')
    logging.debug(f'Keys of the offline_df_dict: {offln_df_dict.keys()}')

    merged_tracks_dict = {}

    for reactor_n, cdf in cdf_dict.items():
        # Debug the problem when reactor is specified as condition in the offline file
        logging.debug(f'Reactor number in the cdf: {reactor_n}')

        # Workaround: in case "conditions" in the offln_df are expressed with letters instead of numbers
        reactor_letters = {1: "A", 2: "B", 3: "C", 4: "D", 5: "E", 6: "F", 7: "G", 8: "H"}

        try:
            offln_df = offln_df_dict[reactor_n]

        except KeyError:
            reactor_l = reactor_letters[reactor_n]
            offln_df = offln_df_dict[reactor_l]

        merged_tracks_dict[reactor_n] = merge_cdf_to_offln_df(cdf, offln_df)

    logging.debug(f'Finished running *MERGE.merge_df_dicts* \n\n')
    return merged_tracks_dict
