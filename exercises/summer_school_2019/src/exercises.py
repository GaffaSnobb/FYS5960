from matplotlib import lines
import numpy as np
import matplotlib.pyplot as plt

bin_size_1 = 0.200    # [MeV].
bin_size_2 = 0.500    # [MeV].

T = 1.238    # [MeV].
E_0 = 0.448   # [MeV].

a = 6.196  # [MeV^-1].
E_1 = 0.942   # [MeV].
A = 56  # Fe-56.

def level_density(bin_size):
    """
    Read level data from file.  Calculate the level density for input
    bin size.

    Parameters
    ----------
    bin_size : float, int
        The bin size in MeV.

    Returns
    -------
    counts : numpy.ndarray
        The number of levels within +- bin_size for all levels.

    counts/bin_size : numpy.ndarray
        The level density.

    levels : numpy.ndarray
        The raw level data.
    """

    levels = np.loadtxt("levels_56Fe.txt")
    levels /= 1e3
    n_levels = len(levels)
    counts = np.zeros(n_levels)

    for i in range(n_levels):
        for j in range(n_levels):
            if (levels[j] >= levels[i] - bin_size) and (levels[j] <= levels[i] + bin_size):
                counts[i] += 1

    return counts, counts/bin_size, levels


def constant_temperature_model(E_x):
    """
    Using global fit parameters (T, E_0) from von Egidy and Bucurescu.
    """
    return np.exp((E_x - E_0)/T)/T


def back_shifted_fermi_gas_model(E_x):
    sigma_j = 0.0146*A**(5/3)*(1 + np.sqrt(1 + 4*a*(E_x - E_1)))/(2*a)
    res = np.exp(2*np.sqrt(a*(E_x - E_1)))
    res /= 12*np.sqrt(2)*sigma_j*a**(1/4)*(E_x - E_1)**(5/4)
    return res


def exercise_1():
    counts_1, density_1, levels = level_density(bin_size_1)
    counts_2, density_2, levels = level_density(bin_size_2)

    plt.step(levels, density_1, label=f"{bin_size_1=}")
    plt.step(levels, density_2, label=f"{bin_size_2=}")
    plt.xlabel("Energy [MeV]")
    plt.ylabel(r"Level density [keV$^{-1}$]")
    plt.legend()
    plt.show()


def exercise_3():
    counts_1, density_1, levels = level_density(bin_size_1)
    counts_2, density_2, levels = level_density(bin_size_2)
    energy_scope = np.linspace(0, 11.197, 1000)
    ctm = constant_temperature_model(energy_scope)
    bsfgm = back_shifted_fermi_gas_model(energy_scope)

    plt.step(levels, density_1, label=f"{bin_size_1=}")
    plt.step(levels, density_2, label=f"{bin_size_2=}")
    plt.plot(energy_scope, ctm, label="CTM")
    plt.plot(energy_scope, bsfgm, label="BSFGM")
    plt.vlines((25/16)/a + E_1, ymin=1e-1, ymax=1e3, linestyles="dashed", color="black")
    plt.title(r"$^{56}$Fe")
    plt.yscale("log")
    plt.xlabel("Energy [MeV]")
    plt.ylabel(r"Level density [MeV$^{-1}$]")
    plt.legend(loc="lower right")
    plt.show()

if __name__ == "__main__":
    # exercise_1()
    exercise_3()