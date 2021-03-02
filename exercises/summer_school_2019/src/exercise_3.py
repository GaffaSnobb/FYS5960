import numpy as np
import matplotlib.pyplot as plt
from exercise_1_levels import level_density
from exercise_2 import ct

# a = 6.196/1000  # [KeV**(-1)].
a = 6.196  # [MeV**(-1)].
E1 = 0.942    # [MeV].
A = 56

def sigma(Ex):
    return 0.0146*A**(5/3)*(1 + np.sqrt(4*a*(Ex - E1)))/(2*a)

def density(Ex):
    return np.exp(2*np.sqrt(a*(Ex - E1)))/(12*np.sqrt(2)*sigma(Ex)*a**(1/4)*(Ex - E1)**(5/4))

if __name__ == "__main__":
    bin_size_1 = 200    # [keV].
    bin_size_2 = 500    # [keV].
    counts_1, density_1, levels = level_density(bin_size_1)
    counts_2, density_2, levels = level_density(bin_size_2)
    levels /= 1000

    plt.step(levels, density_1*1000, label="Discrete level density.")
    plt.plot(levels, density(levels), label="BSFG")
    plt.plot(levels, ct(levels), label="CT")
    plt.yscale("log")
    plt.xlabel("E [MeV]")
    plt.ylabel(r"Density [MeV$^{-1}]$")
    plt.legend()
    plt.show()

