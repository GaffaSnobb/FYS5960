import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def kshell():
    kshell_filename = "kshell_data/summary_V50_gxpf1a.txt"
    bin_width = 0.2
    Ex_max = 100
    Ex_min = 0
    res = ksutil.loadtxt(path=kshell_filename, load_and_save_to_file=True)[0]

    # ksutil.level_plot(
    #     levels = res.levels,
    #     max_spin_states = 10
    # )

    gsf = ksutil.strength_function_average(
        levels = res.levels,
        transitions = res.transitions_BM1,
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = "M1"
    )

    n_bins = int(np.ceil(Ex_max/bin_width))  # NOTE: Why ceil and not floor here?
    Ex_max_adjusted = bin_width*n_bins  # NOTE: Prob. to adjust after the round-off of np.ceil.
    bins = np.linspace(0, Ex_max_adjusted, n_bins + 1)
    bins_middle = (bins[: -1] + bins[1:])/2
    bin_slice = bins_middle[:len(gsf)]
    print(f"{len(bins_middle)=}")
    print(f"{len(gsf)=}")
    print(f"{Ex_max=}")
    print(f"{Ex_max_adjusted=}")
    print(f"{Ex_max/bin_width=}")
    
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
    ax.set_title(r"$^{50}$V")
    plt.show()

if __name__ == "__main__":
    # experimental()
    kshell()