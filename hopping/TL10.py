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
    __df = json_to_dataframe("TL10.json")
    __df['marker'] = ['->', '-<', '-o', '-s', '-^', '-v', '-D', '-X']
    __df['color'] = ['C6', 'C7', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    __df['label'] = ['1', '0.5', r'$6.6\times10^{-2}$', r'$7.4\times10^{-5}$',\
                    r'$8.2\times10^{-6}$', r'$7.3\times10^{-7}$', r'$3.3\times10^{-7}$', '0']
    __eps_R, __eps_A = 127, 45  # permittivity of rutile/anatase
    __m_R, __m_A = 20, 1  # effective mass of rutile/anatase
    __aB_R, __aB_A = 0.34, 2.38  # Bohr radius of rutile/anatase

    @staticmethod
    def get_df():
        """Get full dataframe."""
        return TL10.__df

    @staticmethod
    def get(par=None, sample=None):
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
    def get_samples(start=None, stop=None):
        """Get list of samples."""
        return list(TL10.__df[start:stop].index)

    @staticmethod
    def get_eps(sample=None):
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
    def get_meff(sample=None):
        """Get permittivity.
        :param sample: sample id
        :return: effective mass (parallel, series)
        """
        vR = TL10.get(par='vR', sample=sample)
        if np.isnan(vR):
            meff = {'TL10_5': 7, 'TL10_7.5': 6}[sample]
            return meff, meff  # same value for paralllel and series connection
        return vR*TL10.__m_R+(1-vR)*TL10.__m_A, 1/(vR/TL10.__m_R+(1-vR)/TL10.__m_A)

    @staticmethod
    def read_conduction_data(sample=None) -> 'pd.dataframe':
        # TODO: filenames from dictionary 
        return pd.read_csv(f'datafiles/{sample}.csv', 
            comment='#', index_col=0, header=None, names=['temperature', 'conductance'])

if __name__ == "__main__":
    # print(TL10.get_df())
    df = TL10.read_conduction_data('TL10_10')
    print(df.index.to_numpy())
    # df.plot()
    # plt.show()

