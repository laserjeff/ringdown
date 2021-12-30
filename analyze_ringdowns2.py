"""
      EXCELITAS TECHNOLOGIES
Enabling the future through light.
Created by: Jeff Crawford (7139), jeff.crawford@excelitas.com
Date: 12/30/2021
Department: Metrology
About: This script takes raw ringdown data from Picoscope and fits the data
       to a decay function. The ringdown time is then calculated by multiplying TAU
       by ln(90/10) to get the 90-10 ringdown time in microseconds.
"""
# Non Standard Library:
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
# Standard Library:
from typing import Tuple
import glob
import time
from dataclasses import dataclass


# Decay Function
def decay(x,a,c, max_v) -> float:
    "Decay function where tau (c) * ln(90/10) gives the ringdown time."
    # c = tau
    return a + max_v*np.exp(-x/c) # A + V0*e^(-x/c)

def from_raw_df(data: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    "Given a dataframe created from CSV, returns time normalized dataframe as well as useful parameters labeled in a dictionary."
    # Calculated Ringdown indexes
    start_full = (data['(V)'] == max(data['(V)'])).idxmax()
    end_full = (data['(us)'] >= parameters.time_end).idxmax()
    max_v = max(data['(V)'][start_full:end_full])

    # Normalize ringdown start time to 0
    offset = [data['(us)'][start_full], 0] #[time offset, voltage offset]: Time offset = time where V = max(V), Voltage offset = 0
    data = data.sub(offset) 
    
    # Raw Ringdown indexes
    vrange = max(data['(V)'][start_full:end_full]) - min(data['(V)'][start_full:end_full])
    vmax = max(data['(V)'][start_full:end_full]) - (0.1 * vrange)
    vmin = min(data['(V)'][start_full:end_full]) + (0.1 * vrange)

    start = (data['(V)'][start_full:end_full] <= vmax).idxmax()
    end = (data['(V)'][start_full:end_full] <= vmin).idxmax()
    
    return (data, {'Start': start_full, 'End': end_full, '90 Start': start, '10 End': end, 'Max Voltage': max_v})


def get_calc_ringdown(data: pd.DataFrame, indexes: dict) -> dict:
    '''Fits data to decay function, then calculates the ringdown value
    from tau by multiplying by ln(90/10). Returns ringdown time and standard error (1 stdev).'''

    time = data['(us)'][indexes['Start']:indexes['End']]
    voltage = data['(V)'][indexes['Start']:indexes['End']]
    max_v = indexes['Max Voltage']

    popt, pcov = curve_fit(lambda x, a, c: decay(x,a,c, max_v), 
                           time,
                           voltage,
                           bounds= (parameters.lower_bounds, parameters.upper_bounds)
                           )

    ringdown = popt[1] * np.log(9) #popt[1] is TAU, popt[0] is Voltage offset
    return {'Calc Ringdown': ringdown, 'Error': np.sqrt(np.diag(pcov))[1]}

def get_raw_ringdown(data, indexes) -> pd.DataFrame:
    '''Calculates RAW ringdown through simple subtraction between time at 90% and 10% V values.'''
    raw = data.iloc[indexes['90 Start']:indexes['10 End']]
    raw_ringdown = max(raw['(us)']) - min(raw['(us)'])
    return raw_ringdown

def compile_data(df: pd.DataFrame) -> dict:
    '''Sends DataFrame from CSV data to appropriate functions and returns a row of data containing:
        - Calculated Ringdown time
        - Standard Error of Calculated Ringdown time
        - Raw Ringdown time
        - Difference between calculated and raw ringdown time.'''

    data, indexes = from_raw_df(df)
    calculated_ringdown = get_calc_ringdown(data, indexes)
    r_ringdown = {'Raw Ringdown':get_raw_ringdown(data, indexes)}
    difference = {'Difference': calculated_ringdown['Calc Ringdown'] - r_ringdown['Raw Ringdown']}
    row = {**calculated_ringdown, **r_ringdown, **difference}
    return row


@dataclass
class parameters:
    "Container for parameters"
    data_path: str = 'P:\\JeffC\\ringdown raw data\\11-05-2021 633\\20211105-0002\\'
    time_end: int = 350
    lower_bounds = [-1, 0]
    upper_bounds = [0, 50]

def main():
    '''Glob csv files from data_path and feed into pandas DataFrames,
    then feed DataFrame data into curve_fit to extract TAU.
    Compiles a DataFrame containing the ringdown time from every measurement.
    Returns the MEAN of ringdown time, error, and difference.'''

    start_time = time.time() # Timing how long script takes to run

    data_path = parameters.data_path
    data_files = glob.glob(data_path + '*.csv')

    df_generator = (pd.read_csv(csv, skiprows=1) for csv in data_files)

    data = pd.DataFrame(columns =['Calc Ringdown', 'Error', 'Raw Ringdown', 'Difference'])
    for df in df_generator:
        data = data.append(compile_data(df), ignore_index = True)

    print(data.mean())
    print("")
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()