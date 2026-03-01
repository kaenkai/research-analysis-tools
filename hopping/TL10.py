"""
@author: Karol Kulinowski
Class containing TL10 series data
"""

from lib import json_to_dataframe, KB, fitVRH, fitNNH
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class TL10:
    """Class for handling database of TL10 samples"""
    __df = json_to_dataframe("database/TL10.json")

    # Plot parameters
    __plot_params = {
        'TL10_5': {'m': '->', 'c': 'C6', 'label': '1'},
        'TL10_7.5': {'m': '-<', 'c': 'C7', 'label': '0.5'},
        'TL10_10': {'m': '-o', 'c': 'C0', 'label': r'$6.6\times10^{-2}$'},
        'TL10_15': {'m': '-s', 'c': 'C1', 'label': r'$7.4\times10^{-5}$'},
        'TL10_17.5': {'m': '-^', 'c': 'C2', 'label': r'$8.2\times10^{-6}$'},
        'TL10_20': {'m': '-v', 'c': 'C3', 'label': r'$7.3\times10^{-7}$'},
        'TL10_25': {'m': '-D', 'c': 'C4', 'label': r'$3.3\times10^{-7}$'},
        'TL10_30': {'m': '-X', 'c': 'C5', 'label': '0'}
    }

    # Rutile/anatase properties
    __eps_R, __eps_A = 127, 45  # permittivity of rutile/anatase
    __m_R, __m_A = 20, 1  # effective mass of rutile/anatase
    __aB_R, __aB_A = 0.34, 2.38  # Bohr radius of rutile/anatase

    # Temperature ranges
    __high_temp_ranges = {
        'TL10_5': [180, 210], 'TL10_7.5': [210, 270], 'TL10_10': [220, 270], 'TL10_15': [210, 270],
        'TL10_17.5': [170, 270], 'TL10_20': [220, 280], 'TL10_25': [230, 290], 'TL10_30': [220, 280],
    }
    __low_temp_ranges = {
        'TL10_5': [30, 42], 'TL10_7.5': [30, 40], 'TL10_10': [30, 40], 'TL10_15': [44, 60],
        'TL10_17.5': [100, 120], 'TL10_20': [130, 140], 'TL10_25': [150, 170], 'TL10_30': [170,180],
    }

    @staticmethod
    def plot_params(sample):
        return TL10.__plot_params[sample]

    @staticmethod
    def ht_range(sample):
        return TL10.__high_temp_ranges[sample]

    @staticmethod
    def lt_range(sample):
        return TL10.__low_temp_ranges[sample]

    @staticmethod
    def get_df():
        """Get full dataframe."""
        return TL10.__df

    @staticmethod
    def get(par, sample):
        """Get parameter value for a single sample.
        :return: dataframe cell for a given parameter and sample
        :raise KeyError: if sample or parameter is not found
        """
        try:
            return TL10.__df.at[sample, par]
        except KeyError as err:
            raise KeyError(f'{err}: par={par}, sample={sample}')

    @staticmethod
    def get_range(par=None, start=None, stop=None):
        """Get parameter values for a range of samples.
        :param par: parameter name
        :param start: start index of sample range
        :param stop: stop index of sample range
        :return: parameter values for a given parameter and sample range
        :raise KeyError: if parameter is not found
        """
        try:
            return np.array(TL10.__df[par][start:stop])
        except KeyError as err:
            raise KeyError(f'{err}, wrong parameter name (got {par})')

    @staticmethod
    def get_parameters():
        """Get list of available parameters."""
        return list(TL10.__df.columns)

    @staticmethod
    def get_samples(start, stop):
        """Get list of samples."""
        return list(TL10.__df[start:stop].index)

    @staticmethod
    def get_eps(sample):
        """Get permittivity.
        :param sample: sample id
        :return: permittivity (parallel, series)
        """
        vR = TL10.get(par='vR', sample=sample)
        if np.isnan(vR):
            eps = {'TL10_5': 1, 'TL10_7.5': 5}[sample]
            return eps, eps  # same value for paralllel and series connection
        return vR*TL10.__eps_R+(1-vR)*TL10.__eps_A, 1/(vR/TL10.__eps_R+(1-vR)/TL10.__eps_A)

    @staticmethod
    def get_meff(sample):
        """Get effective mass.
        :param sample: sample id
        :return: effective mass (parallel, series)
        """
        vR = TL10.get(par='vR', sample=sample)
        if np.isnan(vR):
            meff = {'TL10_5': 7, 'TL10_7.5': 6}[sample]
            return meff, meff  # same value for paralllel and series connection
        return vR*TL10.__m_R+(1-vR)*TL10.__m_A, 1/(vR/TL10.__m_R+(1-vR)/TL10.__m_A)

    @staticmethod
    def read_conductivity(sample):
        """Read conductivity data in [1/Ohm/cm]"""
        # TODO: filenames from dictionary
        th = TL10.get(par='thickness_nm', sample=sample)*1e-7
        df = pd.read_csv(f'datafiles/{sample}.csv', 
            comment='#', index_col=0, header=None, names=['temperature', 'conductivity'])/th
        return pd.Series(df['conductivity'], index=df.index)

if __name__ == "__main__":
    # print(TL10.get_df())
    cond = TL10.read_conductivity('TL10_10')
    # print(df.loc[201:99].index.to_numpy())
    print(cond.iat[2])
    # for i, t in enumerate(df.index):
    #     if round(t) == 100: print(df['conductivity'].iat[i])
    # df.plot()
    # plt.show()

