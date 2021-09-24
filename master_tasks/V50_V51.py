import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def kshell_vs_experimental(
    kshell_filename: str,
    experimental_filename: str
    ):
    bin_width = 0.2
    Ex_max = 20
    Ex_min = 0
    res = ksutil.loadtxt(path=kshell_filename, load_and_save_to_file=True)[0]
    
    # res.level_plot(
    #     max_spin_states = 10
    # )
    fig, ax = plt.subplots()

    gsf, bins = ksutil.strength_function_average(
        levels = res.levels,
        transitions = res.transitions_BM1,
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = "M1",
        initial_or_final = "final"
    )
    ax.plot(bins[:50], gsf[:50], label="final")

    gsf, bins = ksutil.strength_function_average(
        levels = res.levels,
        transitions = res.transitions_BM1,
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = "M1",
        initial_or_final = "initial"
    )
    ax.plot(bins[:50], gsf[:50], label="initial")
    
    
    N, Ex, gsf_experimental, gsf_std = np.loadtxt(experimental_filename, skiprows=1, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="Experimental")
    
    ax.set_xlabel("Ex [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.set_title(f"{experimental_filename.split('/')[0]}")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    V50_kshell_filename = "V50/kshell_data/summary_V50_gxpf1a.txt"
    V50_experimental_filename = "V50/V50_gsf.txt"

    V51_kshell_filename = "V51/kshell_data/summary_V51_gxpf1a.txt"
    V51_experimental_filename = "V51/V51_gsf.txt"
    
    kshell_vs_experimental(V50_kshell_filename, V50_experimental_filename)
    # kshell_vs_experimental(V51_kshell_filename, V51_experimental_filename)