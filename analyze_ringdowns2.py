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
def decay(x,a,c, max_v):
    "Decay function where tau (c) * ln(90/10) gives the ringdown time."
    # c = tau
    return a + max_v*np.exp(-x/c) # A + V0*e^(-x/c)

def from_raw_df(data: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    "Given dataframe created from CSV, returns the dataframe with the offset applied as well as key values in a dictionary"

    start_full = (data['(V)'] == max(data['(V)'])).idxmax() #index of first instance of V = max(V) = True
    end_full = (data['(us)'] >= parameters.time_end).idxmax() #index of first instance of time >= user-defined end marker = True
    offset = [data['(us)'][start_full], 0] #offset = [time offset, voltage offset], Voltage offset = 0
    data = data.sub(offset) #Subtract offset from dataframe; [time - time offset, Voltage - 0]
    
    # Finding Range of Voltage values to determine the 90-10 value range      
    vrange = max(data['(V)'][start_full:end_full]) - min(data['(V)'][start_full:end_full])
    vmax = max(data['(V)'][start_full:end_full]) - (0.1 * vrange)
    vmin = min(data['(V)'][start_full:end_full]) + (0.1 * vrange)

    # Start and End positions based on 90-10 value range
    start = (data['(V)'][start_full:end_full] <= vmax).idxmax()
    end = (data['(V)'][start_full:end_full] <= vmin).idxmax()
    max_v = max(data['(V)'][start_full:end_full])
    return (data, {'Start': start_full, 'End': end_full, '90 Start': start, '10 End': end, 'Max Voltage': max_v})


def get_calc_ringdown(data: pd.DataFrame, indexes: dict) -> dict:
    '''Fits data to decay function, then calculates the ringdown value
    from tau by multiplying by ln(90/10). Returns ringdown time and standard error (1 stdev).'''
    time = data['(us)'][indexes['Start']:indexes['End']]
    voltage = data['(V)'][indexes['Start']:indexes['End']]
    max_v = indexes['Max Voltage']

    popt, pcov = curve_fit(lambda x, a, c: decay(x,a,c, max_v), 
                           time,
                           voltage)

    ringdown = popt[1] * np.log(9)
    return {'Ringdown': ringdown, 'Error': np.sqrt(np.diag(pcov))[1]}

def get_raw_ringdown(data, indexes) -> pd.DataFrame:
    raw = data.iloc[indexes['90 Start']:indexes['10 End']]
    raw_ringdown = max(raw['(us)'])-min(raw['(us)'])
    return raw_ringdown

def compile_data(df) -> pd.DataFrame:
    data, indexes = from_raw_df(df)
    calculated_ringdown = get_calc_ringdown(data, indexes)
    r_ringdown = {'Raw Ringdown':get_raw_ringdown(data, indexes)}
    difference = {'Difference': calculated_ringdown['Ringdown']-r_ringdown['Raw Ringdown']}
    rows = {**calculated_ringdown, **r_ringdown, **difference}
    return rows


@dataclass
class parameters:
    "Container for parameters"
    data_path: str = 'P:\\JeffC\\ringdown raw data\\11-05-2021 633\\20211105-0002\\'
    time_end: int = 350

def main():

    start_time = time.time()

    data_path = parameters.data_path
    data_files = glob.glob(data_path + '*.csv')

    df_generator = (pd.read_csv(csv, skiprows=1) for csv in data_files)

    data = pd.DataFrame(columns =['Ringdown', 'Error', 'Raw Ringdown', 'Difference'])
    for df in df_generator:
        data = data.append(compile_data(df), ignore_index = True)

    print(data.mean())
    print("")
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()