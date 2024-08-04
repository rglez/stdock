# Created by gonzalezroy at 8/3/24
"""
This module deals with the spectrum-related tasks
"""
from os.path import join

import pandas as pd
from bruker.api.topspin import Topspin
from bruker.data.nmr import *

import commons as cmn


def get_integral(proton, start, end):
    """
    Get the absolute integral value of a spectrum region

    Args:
        proton: proton object from TopSpin
        start: start value for integration
        end: end value for integration

    Returns:
        summation: absolute integral value
    """
    res = proton.getSpecDataPoints()
    if res.get(EXCEPTION):
        return -1

    startIndex = proton.getIndexFromPhysical(start, 0)
    endIndex = proton.getIndexFromPhysical(end, 0)

    summation = 0
    data = res.get(DATA_POINTS)
    for i in range(startIndex, endIndex):
        summation = summation + data[i]
    return summation


class Spectrum:
    """
    Deal with spectrum-related tasks
    """

    def __init__(self, pdata_path, regions, out_dir='.'):
        # Pdata path & sub folders
        self.path = cmn.check_path(pdata_path)
        self.regions = regions

        # Out dir
        self.root_id = os.sep.join(self.path.split(os.sep)[-3:-1])
        self.out_dir = join(out_dir, self.root_id)

        # TopSpin interface
        self.top = Topspin()
        self.dp = self.top.getDataProvider()

    def integrate(self):
        """
        Integrates the spectrum in the defined regions

        Returns:
            absolute: dataframe of absolute integrals
        """
        # Get TopSpin proton object
        integrals = cmn.recursive_defaultdict()
        proton = self.dp.getNMRData(self.path)
        if proton is None:
            raise ValueError(
                'Something went wrong :( Check the Top Spin connection.')

        # Integrate each region defined in the config file
        for index in self.regions:
            start = self.regions[index]['upper']
            end = self.regions[index]['lower']

            absolute_integral = get_integral(proton, start, end)
            integrals[self.path][index] = absolute_integral

        # Formatting absolute & relative integrals
        absolute = pd.DataFrame(integrals).astype(int)
        absolute['labels'] = [self.regions[x]['label'] for x in self.regions]
        absolute.set_index(absolute.labels, inplace=True)
        del (absolute['labels'])
        return absolute

# =============================================================================
# Debugging Spectrum
# =============================================================================
# pdata_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/12/pdata/11/'
# regions = config_args['std-regions']
# self = Spectrum(pdata_path, regions)

# , total=len(self.commands)
