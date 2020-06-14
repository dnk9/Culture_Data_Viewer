import pandas as pd
import numpy as np

def molarity(cdf: pd.DataFrame, cmpds_properties: pd.DataFrame, usual_cmpds: list):
    """
    Converts the actual concentration from g/L to mol/L.
    Requires a conversion table as input.
    Calculates also CARBON EQUIVALENT molarity.
    Calculates the molarity and C eq molarity for 'added', 'consumed', 'produced'

    Parameters
    ----------
    usual_cmpds : list
        The list of the cmpds which should have the molarity calculated
    cmpds_properties : pd.DataFrame
        A dataframe containing various information for each cmpd. Ex: molar weight, number of carbon
        atoms, density.
        Each cmpd ia on its row.
    cdf : object
        The culture dataframe which will be analysed
    """
    new_cdf = cdf.copy()

    for c in usual_cmpds:  # Selecting only the columns which refers to cmpds
        # Finding the names of the needed columns: 'added' 'consumed' 'produced'
        c_consumed_col = c + "_consumed [g]"
        c_added_col = c + "_add_total [g]"
        c_produced_col = c + "_produced [g]"
        # Creating the names for the new columns
        c_molarity_column = c + " [M]"  # Giving a name to the new column

        c_consumed_M_col = c + "_consumed [M]"
        c_added_M_col = c + "_add_total [M]"
        c_produced_M_col = c + "_produced [M]"


        c_carb_molarity_column = c + " [M C_eq]"  # Giving a name to the new column

        c_consumed_M_Ceq_col = c + "_consumed [M C_eq]"
        c_added_M_Ceq_col = c + "_add_total [M C_eq]"
        c_produced_M_Ceq_col = c + "_produced [M C_eq]"



        cmpd_molarity = cmpds_properties.loc[c]["M [g/mol]"]
        cmpd_C_atoms = cmpds_properties.loc[c]["C Atoms [#]"]

        # SIMPLE Molarity
        c_molarity_values = cdf[c] / cmpd_molarity
        c_consumed_M_val = cdf[c_consumed_col] / cmpd_molarity
        c_added_M_val = cdf[c_added_col] / cmpd_molarity
        c_produced_M_val = cdf[c_produced_col] / cmpd_molarity

        # Ceq Molarity
        c_carb_molarity_values = c_molarity_values * cmpd_C_atoms
        c_consumed_M_Ceq_val = c_consumed_M_val * cmpd_C_atoms
        c_added_M_Ceq_val = c_added_M_val * cmpd_C_atoms
        c_produced_M_Ceq_val = c_produced_M_val * cmpd_C_atoms

        # Saving the values to the cdf
        new_cdf[c_molarity_column] = c_molarity_values
        new_cdf[c_consumed_M_col] = c_consumed_M_val
        new_cdf[c_added_M_col] = c_added_M_val
        new_cdf[c_produced_M_col] = c_produced_M_val

        new_cdf[c_carb_molarity_column] = c_carb_molarity_values
        new_cdf[c_consumed_M_Ceq_col] = c_consumed_M_Ceq_val
        new_cdf[c_added_M_Ceq_col] = c_added_M_Ceq_val
        new_cdf[c_produced_M_Ceq_col] = c_produced_M_Ceq_val

    return new_cdf


def h_clean_balance(cdf: pd.DataFrame, balance_col: str, balance_zero: int):
    new_cdf = cdf.copy()
    # The balance_zero is NOT an index. It refers to the "name" of the row.
    balzero = new_cdf.loc[balance_zero, balance_col]  # What is the value of the balance, at T0?

    new_bal_col_name = balance_col + "ZEROED"  # Ugly name
    new_bal_col_values = (cdf.loc[:, balance_col] - balzero) * (-1)

    new_cdf.loc[:, new_bal_col_name] = new_bal_col_values

    return new_cdf


def mass_added_shot(cdf: pd.DataFrame, cmpd_data_dict) -> pd.DataFrame:
    """
    Calculate the amount of cmpd added to the reactor through shots not given via pump.

    Feed shots are manually injected in the reactor. This information is stored in the time_points xlsx file. This file
    is first directly exported by the program and then re-imported. The file contains the amount of liquid manually
    added (shots) or removed (samples). Before being re-imported it is necessary to specify the kind of liquid added.
    In the tables there is a specific column for that. Each liquid can have more carry more compounds. The 
    cmpd_data_dict translates the liquid name to the contained compounds.

    Parameters
    ----------
    cdf
    cmpd_data_dict

    Returns
    -------

    """

    ignored_liquids = {np.nan, None, "Inoc"}
    added_liquids = set(cdf["Liquid_added"].unique()) - ignored_liquids

    for liquid in added_liquids:
        fltr_r = (cdf.loc[:, "Liquid_added"] == liquid)  # selecting only the rows with the specific liquid addition

        for cmpd in cmpd_data_dict[liquid]:
            cmpd_name = cmpd["cmpd"]
            cmpd_conc = cmpd["Concentration"]

            # Defining new column names
            cmpd_shot_col_name = cmpd_name + "_add_shot [g]"
            cmpd_shot_cumsum_col_name = cmpd_name + "_add_shot_cumsum [g]"

            # Calculate the added mass and save to column
            vol_added = cdf.loc[fltr_r, "Vol_added"] / 1000
            cmpd_shot_mass = vol_added * cmpd_conc

            cdf.loc[fltr_r, cmpd_shot_col_name] = cmpd_shot_mass

            # Calculate the cumulative sum and save to column
            cmpd_shot_mass_filled = cdf[cmpd_shot_col_name].fillna(0)
            cmpd_shot_mass_cumsum = pd.Series.cumsum(cmpd_shot_mass_filled)

            cdf.loc[:, cmpd_shot_cumsum_col_name] = cmpd_shot_mass_cumsum

    return cdf


def mass_added_pump(cdf, cmpd_data_dict):
    # Which are the feeding pumps? Look into the "Feed_pump" column edited in the TP
    fltr_r_feed = cdf["Feed_pump"].notna()
    feeds = cdf[fltr_r_feed]
    try:
        i = 0
        while True:  # For each pump:
            pump = feeds.iloc[i]["Feed_pump"]
            liquid = feeds.iloc[i]["Liquid_added"]
            pump_start_idx = feeds.iloc[i].name

            pump_col_name = f"V{pump}.PV [mL]"
            pump_series = cdf[pump_col_name]

            vol_added = h_series_tare(pump_series, pump_start_idx)  # Tare of the volumes, just for this calculation.

            # For each compound in the liquid fed through the pump #For loop same as shot addition
            for cmpd in cmpd_data_dict[liquid]:
                cmpd_name = cmpd["cmpd"]
                cmpd_conc = cmpd["Concentration"]

                # Defining new column names
                cmpd_feed_pump_col_name = cmpd_name + "_add_feed_pump [g]"
                cmpd_feed_pump_mass = vol_added / 1000 * cmpd_conc

                cdf.loc[:, cmpd_feed_pump_col_name] = cmpd_feed_pump_mass

            i += 1

    except IndexError:
        pass

    return cdf


def mass_added_balance(cdf, cmpd_data_dict):

    fltr_r_feed = cdf["Feed_balance"].notna()
    feeds = cdf[fltr_r_feed]
    try:
        i = 0
        while True:  # For each balance:
            balance = feeds.iloc[i]["Feed_balance"]
            liquid = feeds.iloc[i]["Liquid_added"]
            balance_start_idx = feeds.iloc[i].name

            balance_col_name = f"Bal{balance}.MPV [g]"
            balance_series = cdf[balance_col_name]

            mass_added = h_series_tare(balance_series, balance_start_idx)

            # For each compound in the liquid. #For loop same as shot addition
            for cmpd in cmpd_data_dict[liquid]:
                cmpd_name = cmpd["cmpd"]
                cmpd_conc = cmpd["Concentration"]
                feed_density = cmpd["Density"]

                # Defining new column names
                cmpd_feed_balance_col_name = cmpd_name + "_add_feed_balance [g]"
                cmpd_feed_balance_mass = mass_added / feed_density * cmpd_conc
                cdf.loc[:, cmpd_feed_balance_col_name] = cmpd_feed_balance_mass

            i += 1

    except IndexError:
        pass

    return cdf


def mass_actual(cdf: pd.DataFrame, usual_cmpds: list):
    """
    Calculate the amount of the cmpd in grams.
    Returns the actual mass of the specific cmpd. Multiplies the actual concentration by
    the actual volume.
    """
    new_cdf = cdf.copy()

    for c in usual_cmpds:  # Selecting only the columns which refers to cmpds

        c_mass_column = c + "_actual_mass [g]"  # Name of the new column. Result will be in grams

        c_mass_values = cdf[c] * cdf["Vol_corrected [mL]"] /1000
        new_cdf.loc[:, c_mass_column] = c_mass_values

    return new_cdf


def mass_added_initial(cdf, usual_cmpds):
    cdf = cdf.copy()
    for c in usual_cmpds:
        mass_column = c + '_actual_mass [g]'
        initial_mass_col_name = c + '_add_initial [g]'
        idx = cdf[mass_column].first_valid_index()
        if idx is None:
            initial_value = 0
        else:
            initial_value = cdf.iloc[idx][mass_column]
        cdf[initial_mass_col_name] = initial_value
    return cdf


def mass_added_total(cdf, added_cmpds_dict):
    """Calculate the total amount of each added compound. It is the sum of all the '_add_' columns."""
    cdf = cdf.copy()
    #In this function I can check if the mass added from pump and balance are ok
    for cmpd, col_list in added_cmpds_dict.items():
        col_name_mass_added_total = cmpd + "_add_total [g]"
        # Warning! Here it will sum the mass calculated from the balance to the one calculated from the pump!
        # Note: the columns are interpolated now, but this should be done earlier!
        col_value = cdf.loc[:, col_list].interpolate().sum(axis=1)
        cdf[col_name_mass_added_total] = col_value
    return cdf



def mass_removed_total(cdf, removed_cmpds_dict):
    """Calculate the total amount of each added compound. It is the sum of all the '_add_' columns."""
    cdf = cdf.copy()
    for cmpd, col_list in removed_cmpds_dict.items():
        col_name_mass_removed_total = cmpd + "_rem_total [g]"
        # Note: the columns are interpolated now, but this should be done earlier!
        col_value = cdf.loc[:, col_list].interpolate().sum(axis=1)
        cdf[col_name_mass_removed_total] = col_value
    return cdf


def mass_removed_sampling(cdf, usual_cmpds):
    for i in usual_cmpds:
        # Not always the values from the TP and the ones from the sampling are on
        # the same line. To be sure to avoid errors, here I am interpolating and
        # filling the NA before the first sample and after the last one. The
        # distance between the values from the different columns should be anyway
        # short, and so the interpolation caused deviation.
        i_int = cdf[i].interpolate().fillna(method="bfill").fillna(method="ffill")
        i_mass_removed_sampling = i_int * cdf["Vol_removed"] / 1000
        i_mass_removed_sampling_cumsum = pd.Series.cumsum(i_mass_removed_sampling.fillna(0))
        # Mass removed with each sample
        col_name = i + "_rem_sample [g]"
        cdf[col_name] = i_mass_removed_sampling

        # Total mass removed, cumulative sum
        col_name_cumsum = i + "_rem_sample_cumsum [g]"
        cdf[col_name_cumsum] = i_mass_removed_sampling_cumsum

    return cdf


def mass_consumed_produced(cdf, removed_cmpds_dict):
    """Create the 'consumed' and 'produced' columns for each compound. """
    # Using 'removed_cmpds_dict' because for sure has all the columns needed
    # maybe also 'usual_cmpds'?
    cdf = cdf.copy()
    for k in removed_cmpds_dict.keys():
        add_col = k + '_add_total [g]'
        rem_col = k + '_rem_total [g]'
        act_col = k + '_actual_mass [g]'
        cons_col = k + '_consumed [g]'
        prod_col = k + '_produced [g]'
        cdf[cons_col] = cdf[add_col] - cdf[act_col] - cdf[rem_col]
        cdf[prod_col] = - cdf[add_col] + cdf[act_col] + cdf[rem_col]
    return cdf


def h_series_tare(series, index) -> pd.Series:
    # In case there is no value at the index, take the last value
    zero = series[index - 50:index + 1].fillna(method="ffill")[index]

    series = (series - zero).abs()
    series[:index] = 0
    return series


def h_extract_cmpd_list(columns, cmpd_set) -> set:
    """ Make a list of the columns which refers to valid compounds. """
    return {i for i in columns if i.casefold() in cmpd_set}


def h_list_added_cmpds(columns) -> dict:
    """ Make a dict with the compounds which are fed. Key is the compound name, value is the list of the columns which
    needs to be summed"""
    # "Valid" columns have this format: 'compoundname_add_feed_balance [g]' or 'compoundname_add_shot [g]'
    # TODO: here! focus only on the cumsum (in cases where both istantaneous and cumsum are available)
    # This list/dict is made only to ease the generation of the "added_total" column. "shot" non cumsum is not needed
    # for this calculation.
    # TODO: decide which to take between BALANCE and PUMP (both should give the same value, but are now summed)
    cols_with_valid_name = {i for i in columns if "_add_" in i}
    added_cmpds = {}
    for c in cols_with_valid_name:
        cmpd = c.split("_")[0]
        if cmpd not in added_cmpds.keys():
            added_cmpds[cmpd] = [c]
        else:
            added_cmpds[cmpd].append(c)

    return added_cmpds


def h_list_removed_cmpds(columns) -> dict:
    """ Make a dict with the compounds which are removed. Key is the compound name, value is the list of the columns
    which needs to be summed"""
    # "Valid" columns have this format: 'compoundname_rem_sample [g]' or 'compoundname_rem_pump [g]'
    # TODO: here! focus only on the cumsum (in cases where both istantaneous and cumsum are available)
    # This list/dict is made only to ease the generation of the "removed_total" column. "sample" non cumsum is not
    # needed for this calculation.
    cols_with_valid_name = {i for i in columns if "_rem_" in i}
    removed_cmpds = {}
    for c in cols_with_valid_name:
        cmpd = c.split("_")[0]
        if cmpd not in removed_cmpds.keys():
            removed_cmpds[cmpd] = [c]
        else:
            removed_cmpds[cmpd].append(c)

    return removed_cmpds


def yields(cdf, usual_cmpds, substrates):
    cdf = cdf.copy()
    for cmpd in usual_cmpds:
        if cmpd in substrates:
            pass
        else:
            cmpd_produced_col = cmpd + '_produced [g]'
            cmpd_produced_val = cdf[cmpd_produced_col]

            cmpd_produced_M_col = cmpd + '_produced [M]'
            cmpd_produced_M_val = cdf[cmpd_produced_M_col]

            cmpd_produced_M_C_eq_col = cmpd + '_produced [M C_eq]'
            cmpd_produced_M_C_eq_val = cdf[cmpd_produced_M_C_eq_col]

            for sub in substrates:
                # Yield in g/g
                sub_consumed_col = sub + '_consumed [g]'
                sub_total_col = sub + '_add_total [g]'
                yield_consumed_col = cmpd + '_yield_consumed_' + sub + ' [g/g]'
                yield_total_col = cmpd + '_yield_total_' + sub + ' [g/g]'
                #print(f"yield g/g consumed = {cmpd_produced_col}/{sub_consumed_col} ")
                #print(f"yield g/g total = {cmpd_produced_col}/{sub_total_col} ")
                sub_consumed_val = cdf[sub_consumed_col]
                sub_total_val = cdf[sub_total_col]

                yield_consumed_val = cmpd_produced_val / sub_consumed_val
                yield_total_val = cmpd_produced_val / sub_total_val

                cdf[yield_consumed_col] = yield_consumed_val
                cdf[yield_total_col] = yield_total_val

                # Yield in M/M
                sub_consumed_M_col = sub + '_consumed [M]'
                sub_total_M_col = sub + '_add_total [M]'
                yield_consumed_M_col = cmpd + '_yield_consumed_' + sub + ' [M/M]'
                yield_total_M_col = cmpd + '_yield_total_' + sub + ' [M/M]'
                #print(f"yield M/M consumed = {cmpd_produced_col}/{sub_consumed_col} ")
                #print(f"yield M/M total = {cmpd_produced_col}/{sub_total_col} ")
                sub_consumed_M_val = cdf[sub_consumed_M_col]
                sub_total_M_val = cdf[sub_total_M_col]

                yield_consumed_M_val = cmpd_produced_M_val / sub_consumed_M_val
                yield_total_M_val = cmpd_produced_M_val / sub_total_M_val

                cdf[yield_consumed_M_col] = yield_consumed_M_val
                cdf[yield_total_M_col] = yield_total_M_val

                # Yield in M Ceq/ M Ceq
                sub_consumed_M_C_eq_col = sub + '_consumed [M C_eq]'
                sub_total_M_C_eq_col = sub + '_add_total [M C_eq]'
                yield_consumed_M_C_eq_col = cmpd + '_yield_consumed_' + sub + ' [M/M (C_eq)]'
                yield_total_M_C_eq_col = cmpd + '_yield_total_' + sub + ' [M/M (C_eq)]'

                sub_consumed_M_C_eq_val = cdf[sub_consumed_M_C_eq_col]
                sub_total_M_C_eq_val = cdf[sub_total_M_C_eq_col]

                yield_consumed_M_C_eq_val = cmpd_produced_M_C_eq_val / sub_consumed_M_C_eq_val
                yield_total_M_C_eq_val = cmpd_produced_M_C_eq_val / sub_total_M_C_eq_val

                cdf[yield_consumed_M_C_eq_col] = yield_consumed_M_C_eq_val
                cdf[yield_total_M_C_eq_col] = yield_total_M_C_eq_val

    return cdf


def rates(cdf, cmpds_substrates: set, cmpds_measured: set):

    cdf = cdf.copy()
    fltr_r_samples = (cdf["Sample"] == True)
    smp = cdf[fltr_r_samples]

    # Making a list of the new columns which I am calculating. Will be necessary later when I merge the "sample" df
    # back to the original
    new_cols_names = []

    for cmpd in cmpds_measured:

        if cmpd in cmpds_substrates:
            status = '_consumed '
            status_c = '_consumption '
        else:
            status = '_produced '
            status_c = '_production '
            # If cmpd is not defined as a substrate, assume it is a product

        rates_avg_gl = [np.nan]
        rates_inst_gl = [np.nan]
        rates_avg_M = [np.nan]
        rates_inst_M = [np.nan]

        t = 1
        while t < smp.shape[0]:

            tdelta_avg = smp.iloc[t]['InoculationTime []'] / pd.Timedelta(hours=1)
            tdelta_inst = (smp.iloc[t]['InoculationTime []'] - smp.iloc[t - 1]['InoculationTime []']) / pd.Timedelta(
                hours=1)

            vol_avg = smp.iloc[0:t + 1]["Vol_corrected [mL]"].mean()
            vol_inst = smp.iloc[t - 1:t + 1]["Vol_corrected [mL]"].mean()

            cmpd_avg_g = smp.iloc[t][f'{cmpd}{status}[g]']
            cmpd_avg_M = smp.iloc[t][f'{cmpd}{status}[M]']

            cmpd_inst_g = (smp.iloc[t][f'{cmpd}{status}[g]'] - smp.iloc[t - 1][f'{cmpd}{status}[g]'])
            cmpd_inst_M = (smp.iloc[t][f'{cmpd}{status}[M]'] - smp.iloc[t - 1][f'{cmpd}{status}[M]'])

            # Average rate, g/L/h
            rate_avg_gl = (cmpd_avg_g / tdelta_avg) / (vol_avg / 1000)
            rates_avg_gl.append(rate_avg_gl.round(2))

            # Instant rate, g/L/h
            rate_inst_gl = (cmpd_inst_g / tdelta_inst) / (vol_inst / 1000)
            rates_inst_gl.append(rate_inst_gl.round(2))

            # Average rate, M/h
            rate_avg_M = (cmpd_avg_M / tdelta_avg) / (vol_avg / 1000)
            rates_avg_M.append(rate_avg_M.round(6))

            # Instant rate, M/h
            rate_inst_M = (cmpd_inst_M / tdelta_inst) / (vol_inst / 1000)
            rates_inst_M.append(rate_inst_M.round(6))

            t += 1

        rates_avg_gl_col_name = f'{cmpd}{status_c} rate (avg) [g/L/h]'
        rates_inst_gl_col_name = f'{cmpd}{status_c} rate (inst) [g/L/h]'
        rates_avg_M_col_name = f'{cmpd}{status_c} rate (avg) [M/h]'
        rates_inst_M_col_name = f'{cmpd}{status_c} rate (inst) [M/h]'

        new_cols_names.extend([rates_avg_gl_col_name, rates_inst_gl_col_name, rates_avg_M_col_name, rates_inst_M_col_name])

        smp[rates_avg_gl_col_name] = rates_avg_gl
        smp[rates_inst_gl_col_name] = rates_inst_gl
        smp[rates_avg_M_col_name] = rates_avg_M
        smp[rates_inst_M_col_name] = rates_inst_M

    merged_cdf = cdf.merge(smp[new_cols_names], how = "left", left_index=True, right_index=True)

    return merged_cdf


def selectivity(cdf, cmpds_products):
    cdf = cdf.copy()

    for c in cmpds_products:
        # How many grams of c?
        c_prod_g_col_name = f'{c}_produced [g]'
        c_prod_M_col_name = f'{c}_produced [M]'
        c_prod_M_ceq_col_name = f'{c}_produced [M C_eq]'
        for d in cmpds_products - {c}:
            # For how many grams of sub?
            d_prod_g_col_name = f'{d}_produced [g]'
            d_prod_M_col_name = f'{d}_produced [M]'
            d_prod_M_ceq_col_name = f'{d}_produced [M C_eq]'

            sel_g_col_name = f'{c}_selectivity_{d} [g/g]'
            sel_M_col_name = f'{c}_selectivity_{d} [M/M]'
            sel_M_ceq_col_name = f'{c}_selectivity_{d} (C eq)[M/M]'

            # values
            sel_g_col_val = cdf[c_prod_g_col_name] / cdf[d_prod_g_col_name]
            sel_M_col_val = cdf[c_prod_M_col_name] / cdf[d_prod_M_col_name]
            sel_M_ceq_col_val = cdf[c_prod_M_ceq_col_name] / cdf[d_prod_M_ceq_col_name]

            # apply values
            cdf[sel_g_col_name] = sel_g_col_val
            cdf[sel_M_col_name] = sel_M_col_val
            cdf[sel_M_ceq_col_name] = sel_M_ceq_col_val

    return cdf


def co2calc(cdf, choice="PV"):
    cdf = cdf.copy()  # Calculating CO2 grams

    # choice = 'PV' #or 'Out' or 'SV'
    outgas_CO2_prcnt_col_name = 'XCO2.Out [%]'
    outgas_CO2_prcnt_col_val = cdf[outgas_CO2_prcnt_col_name]

    # To fill the gaps and smooth the curve
    outgas_CO2_prcnt_smoothed_val = outgas_CO2_prcnt_col_val.interpolate().rolling(window=200, center=True).mean()
    outgas_flow_col_name = f'F.{choice} [sL/h]'
    outgas_flow_col_val = cdf[outgas_flow_col_name]
    # Smoothing the curves
    outgas_flow_smoothed = outgas_flow_col_val.interpolate().rolling(window=200, center=True).mean()

    # outgas_CO2_flow_val = outgas_CO2_prcnt_col_val * outgas_flow_col_val /100
    outgas_CO2_flow_val = outgas_CO2_prcnt_smoothed_val * outgas_flow_smoothed / 100

    # - in
    ingas_CO2_flow_col_name = 'FCO2.PV [sL/h]'
    ingas_CO2_flow_val = cdf[ingas_CO2_flow_col_name]
    ingas_CO2_flow_smoothed = ingas_CO2_flow_val.interpolate().rolling(window=200, center=True).mean()

    # Difference -> now I have sL/h of CO2, time to go to g/h
    produced_CO2_flow_val = outgas_CO2_flow_val - ingas_CO2_flow_smoothed

    # Moles, grams produced /h
    produced_CO2_Mh_val = produced_CO2_flow_val / 22.414
    produced_CO2_gh_val = produced_CO2_Mh_val * 44.009

    # Calculate total CO2 grams produced: cumulative sum of all the measures.
    # As measurments are taken every 30 seconds, there will be 120 measurements in one hour.
    # As the measurements refer to h^-1, I have to divide the cumsum by 120
    # TODO: this is very reactor dependent -> calculate automatically the value

    CO2_g_prod_cumsum = pd.Series.cumsum(produced_CO2_gh_val) / 120
    CO2_M_prod_cumsum = pd.Series.cumsum(produced_CO2_Mh_val) / 120


    # Saving to cdf
    cdf['CO2_L_prod_rate [sL/h]'] = produced_CO2_flow_val

    cdf['CO2_g_prod_rate [g/h]'] = produced_CO2_gh_val
    cdf['CO2_M_prod_rate [M/h]'] = produced_CO2_Mh_val
    cdf['CO2_g_prod_cumsum [g]'] = CO2_g_prod_cumsum
    cdf['CO2_M_prod_cumsum [M]'] = CO2_M_prod_cumsum

    return cdf


def complete_analysis(cdf, cmpds_substrates, cmpds_properties, cmpd_set, cmpd_data_dict):

    a = cdf.copy()
    cols = a.columns

    # cmpds_substrates = {'Glycerol'}
    cmpds_measured = h_extract_cmpd_list(cols, cmpd_set)
    cmpds_products = cmpds_measured - cmpds_substrates

    b = mass_actual(a, cmpds_measured)
    c = mass_added_initial(b, cmpds_measured)
    # d = mass_added_shot(c, cmpd_data_dict) # Will give problems if no shots
    r = mass_added_pump(c, cmpd_data_dict)
    t = mass_added_balance(r, cmpd_data_dict)

    cmpds_added = h_list_added_cmpds(t.columns)

    u = mass_added_total(t, cmpds_added)
    v = mass_removed_sampling(u, cmpds_measured)

    cmpds_removed = h_list_removed_cmpds(v.columns)

    z = mass_removed_total(v, cmpds_removed)

    aa = mass_consumed_produced(z, cmpds_removed)
    ab = molarity(aa, cmpds_properties, cmpds_measured)
    ac = yields(ab, cmpds_measured, cmpds_substrates)
    ad = rates(ac, cmpds_substrates, cmpds_measured)
    ae = selectivity(ad, cmpds_products)
    af = co2calc(ae)

    return af
