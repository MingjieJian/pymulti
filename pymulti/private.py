import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import shutil
import re
import pkg_resources
import subprocess
from scipy.signal import savgol_filter
try:
    from mendeleev import element
except ImportError:
    print('mendeleev is not imported')
    
package_path = pkg_resources.resource_filename('pymulti', '')
    
def calculate_delta(z):
    '''
    A function to calculate delta of an array.

    Parameters
    ----------
    z : np.array
        The geometrical depth of atmosphere grid.

    Returns
    -------
    dz : np.array
        An array of delta z, with the first element as 0.
    '''

    if type(z) != 'numpy.ndarray':
        z = np.array(z)
    dz = [0]
    for i in range(len(z)-1):
        dz.append(z[i+1] - z[i])
    return np.array(dz)

def load_asplund_abundance():
    return pd.read_csv(package_path + '/files/asplund09.csv')

asplund09 = load_asplund_abundance()