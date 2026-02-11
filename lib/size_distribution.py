import scipy
import numpy as np
from abc import ABC


class SizeDistribution(ABC):
    def __init__(self, distribution_function, r_min, r_max):
        self._pdf = distribution_function
        self.r_min = r_min
        self.r_max = r_max
        self._normalisation = scipy.integrate.quad(self._pdf, r_min, r_max)[0]
    
    def pdf(self, r):
        return np.where(
            r < self.r_min, 0,
            np.where(r > self.r_max, 0,
                     self._pdf(r) / self._normalisation))

    def nth_moment(self, n, r_min=None, r_max=None):
        r_min = max(r_min or self.r_min, self.r_min)
        r_max = min(r_max or self.r_max, self.r_max)
        return scipy.integrate.quad(lambda x: x**n * self.pdf(x), r_min, r_max)[0]


class Revell2021GammaSizeDistribution(SizeDistribution):
    """A gamma distribution with shape parameter 2 and scale parameter 15 μm,
    derived from various observations in Revell et al. (2021)."""
    def __init__(self, r_min, r_max):
        super().__init__(lambda r: scipy.stats.gamma.pdf(2*r, 2, scale=15), r_min, r_max)


class Leusch2023PowerLawSizeDistribution(SizeDistribution):
    """A power law distribution with exponent -1.52, derived from the data of
    Leusch et al. (2023) for outdoor airborne microplastics."""
    def __init__(self, r_min, r_max):
        super().__init__(lambda r: (2*r) ** (-1.52), r_min, r_max)


class PowerLawSizeDistribution(SizeDistribution):
    """A power law distribution with a custom exponent."""
    def __init__(self, exponent, r_min, r_max):
        super().__init__(lambda r: (2*r) ** exponent, r_min, r_max)