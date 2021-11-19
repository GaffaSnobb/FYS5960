import time
import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def compare_V50_experimental():
    fig, ax = plt.subplots()
    V50_experimental_filename = "V50_gsf.txt"
    N, Ex, gsf_experimental, gsf_std = np.loadtxt(V50_experimental_filename, skiprows=1, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="A.C. Larsen et al., Phys. Rev. C 73, 064301 (2006)")

    bin_width = 0.2
    Ex_min = 5
    Ex_max = 20
    multipole_type = "M1"
    
    res = ksutil.loadtxt(
        path = "gxpf1a/400_states_corrected/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label="GXPF1A M1 400 states")

    ax.set_ylim([2.0610052684616102e-09, 1.6590428320549444e-08])
    ax.set_xlabel(r"E$_{\gamma}$ [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.set_title(r"$^{50}$V")
    ax.legend(fontsize=10)
    fig.savefig("v50_gsf_kshell_vs_experimental_gxpf1a_400_states.png", dpi=300)
    plt.show()

def compare_100_200_400_states():
    """
    Compare the gamma strength function using 100, 200 and 400 of the
    lowest energy levels for all spins from 0 to 14. Compare with
    experimental data for V50.
    """
    fig, ax = plt.subplots()
    V50_experimental_filename = "V50_gsf.txt"
    N, Ex, gsf_experimental, gsf_std = np.loadtxt(V50_experimental_filename, skiprows=1, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="A.C. Larsen et al., Phys. Rev. C 73, 064301 (2006)")

    bin_width = 0.2
    Ex_min = 5
    Ex_max = 20
    multipole_type = "M1"

    res = ksutil.loadtxt(
        path = "gxpf1a/400_states/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label="GXPF1A 400 states")
    
    res = ksutil.loadtxt(
        path = "gxpf1a/200_states/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label="GXPF1A 200 states")
    
    res = ksutil.loadtxt(
        path = "gxpf1a/100_states/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label="GXPF1A 100 states")
    
    ax.set_ylim([2.0610052684616102e-09, 1.6590428320549444e-08])
    # ax.set_yscale("log")
    ax.set_xlabel(r"E$_{\gamma}$ [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.set_title(r"$^{50}$V")
    ax.legend(fontsize=10)
    fig.savefig("v50_gsf_kshell_vs_experimental_gxpf1a_100_200_400_states.png", dpi=300)
    plt.show()

def compare_gyroscopic_factors():
    fig, ax = plt.subplots()
    V50_experimental_filename = "V50_gsf.txt"
    N, Ex, gsf_experimental, gsf_std = np.loadtxt(V50_experimental_filename, skiprows=1, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="A.C. Larsen et al., Phys. Rev. C 73, 064301 (2006)")

    bin_width = 0.2
    Ex_min = 5
    Ex_max = 20
    multipole_type = "M1"
    
    res = ksutil.loadtxt(
        path = "gxpf1a/400_states/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    label = r"$g_{s,p}, g_{s,n} = $"
    label += f"{res.parameters['gs'][0]/ksutil.GS_FREE_PROTON:.0f}"
    label += r"$g_{s free, p}, $"
    label += f"{res.parameters['gs'][1]/ksutil.GS_FREE_NEUTRON:.0f}"
    label += r"$g_{s free, n}$"
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label=label)

    res = ksutil.loadtxt(
        path = "gxpf1a/400_states_corrected/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    label = r"$g_{s,p}, g_{s,n} = $"
    label += f"{res.parameters['gs'][0]/ksutil.GS_FREE_PROTON:.1f}"
    label += r"$g_{s free, p}, $"
    label += f"{res.parameters['gs'][1]/ksutil.GS_FREE_NEUTRON:.1f}"
    label += r"$g_{s free, n}$"
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label=label)
    
    ax.set_xlabel(r"E$_{\gamma}$ [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.set_title(r"$^{50}$V")
    ax.legend(fontsize=10)
    fig.savefig("v50_gsf_kshell_vs_experimental_gxpf1a_400_states_gs_comparison.png", dpi=300)
    plt.show()

def E2():
    fig, ax = plt.subplots()
    V50_experimental_filename = "V50_gsf.txt"
    N, Ex, gsf_experimental, gsf_std = np.loadtxt(V50_experimental_filename, skiprows=1, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="A.C. Larsen et al., Phys. Rev. C 73, 064301 (2006)")

    bin_width = 0.2
    Ex_min = 5
    Ex_max = 20
    multipole_type = "M1"
    
    res = ksutil.loadtxt(
        path = "gxpf1a/400_states_corrected/summary_V50_gxpf1a.txt",
        load_and_save_to_file = True
    )[0]
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label=multipole_type)

    multipole_type = "E2"
    
    bins, gsf = res.gsf(
        bin_width = bin_width,
        Ex_min = Ex_min,
        Ex_max = Ex_max,
        multipole_type = multipole_type,
        plot = False
    )
    ax.plot(bins[0:50], gsf[0:50], label=multipole_type)
    
    ax.set_xlabel(r"E$_{\gamma}$ [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.set_title(r"$^{50}$V")
    ax.legend(fontsize=10)
    # fig.savefig("aidz.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    # compare_100_200_400_states()
    # compare_V50_experimental()
    # compare_gyroscopic_factors()
    E2()