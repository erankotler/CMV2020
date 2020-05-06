'''
Concentration of functions for data normalization and manipulations
'''

import numpy as np
import pandas as pd


def standardize_df(df, column='all'):
    ''' 
    Standardize each column of a data frame (change to z scores) & return standardized data frame
    if a specific column name (string) is given - standardize only this column
    '''
    if column =='all':
        std_df = pd.DataFrame(index=df.index)
        for c in df.columns:
            mean = df[c].mean()
            stdev = np.std(df[c])
            std_df[c] = (df[c]-mean)/float(stdev)
    else:
        std_df = df
        mean = df[column].mean()
        stdev = np.std(df[column])
        std_df[column] = (df[column]-mean)/float(stdev)    
        
    return std_df

