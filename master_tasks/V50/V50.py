import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def kshell():
    kshell_filename = "kshell_data/summary_V50_gxpf1a.txt"
    bin_width = 0.2
    E_max = 100
    E_min = 0
    res = ksutil.loadtxt(path=kshell_filename)[0]

    # ksutil.level_plot(
    #     levels = res.levels,
    #     max_spin_states = 10
    # )

    # print(res.transitions[0])

    Egs = res.levels[0, 0]
    res.transitions[:, 2] += Egs
    # res.transitions_BM1[:, 2] += Egs
    # res.transitions_BE2[:, 2] += Egs
    # print(res.transitions[:, 2])
    print(res.levels[:, 0].shape)

    gsf = ksutil.strength_function_average(
        levels = res.levels,
        transitions = res.transitions_BM1,
        Jpi_list = ksutil.create_jpi_list(spins=res.levels[:, 1], parities=res.levels[:, 2]),
        bin_width = bin_width,
        Ex_min = E_min,
        Ex_max = E_max,
        multipole_type = "M1"
    )
    # print(f"{gsf=}")

    n_bins = int(np.ceil(E_max/bin_width))
    E_max_adjusted = bin_width*n_bins
    bins = np.linspace(0, E_max_adjusted, n_bins + 1)
    bins_middle = (bins[0: -1] + bins[1:])/2
    bin_slice = bins_middle[0:len(gsf)]
    
    fig, ax = plt.subplots()
    ax.plot(bin_slice[:50], gsf[:50])
    plt.show()


def experimental():
    experimental_filename = "V50_gsf.txt"
    N, Ex, gsf, gsf_std = np.loadtxt(experimental_filename, skiprows=1, unpack=True)

    fig, ax = plt.subplots()
    ax.errorbar(Ex, gsf, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="Experimental")
    ax.set_xlabel("Ex [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    # experimental()
    kshell()