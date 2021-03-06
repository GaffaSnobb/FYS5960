import sys, os
import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def level_density(levels, bin_size):
    """
    UNUSED
    Calculate the level density for a given bin size.

    Parameters
    ----------
    levels : numpy.ndarray
        Energy levels.

    bin_size : int, float
        The interval of which to calculate the density.

    Returns
    -------
    counts/bin_size : numpy.ndarray
        The level density.
    """
    bins = np.arange(0, levels[-1] + bin_size, bin_size)
    counts = np.zeros(len(bins))

    for i in range(len(bins) - 1):
        counts[i] = np.sum(bins[i] <= levels[levels < bins[i + 1]])

    return (counts/bin_size)[:-1], bins[1:]

def mean_transition_strength(B_M1):
    """
    UNUSED
    Calculate the mean value of the reduced transition probability for
    each energy level.

    Parameters
    ----------
    B_M1 : list
        Nested list of initial energy level and reduced transition
        probability of said level when decaying by an E_gamma energy
        photon.  Structure: [[E_i, b_m1, E_gamma], ...].

    Returns
    -------
    B_M1_mean : numpy.ndarray
        The mean reduced transition probability for each energy level.

    Raises
    ------
    AssertionError
        If MAX_ITERATIONS is reached before the file is read through.
    """
    current_idx = 0         # Index of E_i values in B_M1.
    B_M1_mean = []          # Mean of reduced transition probabilities.
    B_M1_mean_level = []    # Corresponding energy levels.
    MAX_ITERATIONS = int(1e5)
    iteration_counter = 0

    while (current_idx < B_M1.shape[0]) and (iteration_counter < MAX_ITERATIONS):
        """
        Find all elements in the 0th column of B_M1 (E_i values) which
        are identical to the 'current_idx'th value of the same column.
        Then, take the mean of the b_m1 values in the same rows as the
        identical E_i values.

        Update 'current_idx' to to start at the next new E_i value.
        """
        res_idx = np.where(B_M1[:, 0] == B_M1[current_idx, 0])
        B_M1_mean.append(np.mean(B_M1[:, 1][res_idx]))
        B_M1_mean_level.append(B_M1[current_idx, 0])
        
        current_idx = res_idx[0][-1] + 1    # Choose the next new E_i value.
        iteration_counter += 1

    msg = f"Warning! Maximum number of iterations ({MAX_ITERATIONS}) reached!"
    msg += " The data file is either very large or something went wrong."
    success = iteration_counter != MAX_ITERATIONS
    assert success, msg


    return np.array(B_M1_mean), np.array(B_M1_mean_level)

def plot_gsf(isotope_name):
    """
    Plot the gamma strength function for a single isotope.

    isotope_name : string
        Examples: S24, Ne30.

    Raises
    ------
    ValueError
        If isotope_name cannot be found in the calculated data
        files.
    """
    fname = None
    
    for fnames in self.fnames_combined:
        for i in range(len(fnames)):
            if isotope_name in fnames[i][0]:
                fname = fnames[i][0]

    if fname is None:
        msg = f"Isotope name '{isotope_name}' is not a valid name."
        raise ValueError(msg)

    res = ksutil.loadtxt(self.directory + fname)

    _, ax = plt.subplots()

    Jpi_list = ksutil.create_jpi_list(res.levels[:, 1], None)
    E_gs = res.levels[0, 0]
    res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with J??rgen.

    gsf = ksutil.strength_function_average(
        levels = res.levels,
        transitions = res.transitions,
        Jpi_list = Jpi_list,
        bin_width = self.bin_width,
        Ex_min = self.Ex_min,    # [MeV].
        Ex_max = self.Ex_max,    # [MeV].
        multipole_type = "M1"
    )

    bin_slice = self.bins_middle[0:len(gsf)]
    ax.plot(bin_slice, gsf, label=fname)
    ax.legend()
    ax.set_xlabel(r"$E_{\gamma}$ [MeV]")
    ax.set_ylabel(r"gsf [MeV$^{-3}$]")
    plt.show()

def recreate(directory, filter_=None):
    """
    Recreate figure 6 from JEM's paper.
    """
    if filter_ is None:
        elements = [
            "oxygen", "fluorine", "neon", "argon"
        ]
    else:
        elements = [filter_]
    
    bin_width = 0.2
    E_max = 46
    Ex_min = 0
    Ex_max = E_max
    n_bins = int(np.ceil(E_max/bin_width))
    E_max_adjusted = bin_width*n_bins
    bins = np.linspace(0, E_max_adjusted, n_bins + 1)
    bins_middle = (bins[0: -1] + bins[1:])/2
    fig, ax = plt.subplots()

    for element in elements:

        res_list = ksutil.loadtxt(
            path = directory,
            is_directory = True,
            filter_ = element,
            load_and_save_to_file = True
        )

        ratios = []
        for res in res_list:

            Jpi_list = ksutil.create_jpi_list(
                spins = res.levels[:, 1],
                parities = res.levels[:, 2]
            )
            E_gs = res.levels[0, 0]

            try:
                res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with J??rgen.
            except IndexError:
                print(f"File {res.fname_summary} skipped! Too few / no energy levels are present in this data file.")
                ratios.append(None) # Maintain correct list length for plotting.
                continue
            
            gsf = ksutil.strength_function_average(
                levels = res.levels,
                transitions = res.transitions,
                Jpi_list = Jpi_list,
                bin_width = 0.2,
                Ex_min = Ex_min,    # [MeV].
                Ex_max = Ex_max,    # [MeV].
                multipole_type = "M1"
            )

            # Sum gsf for low and high energy range and take the ratio.
            bin_slice = bins_middle[0:len(gsf)]
            low_idx = (bin_slice <= 2)
            high_idx = (bin_slice <= 6) == (2 <= bin_slice)
            low = np.sum(gsf[low_idx])
            high = np.sum(gsf[high_idx])
            low_high_ratio = low/high
            ratios.append(low_high_ratio)

        ax.plot(range(8, 20+1), ratios, "--.", label=element)
        ax.set_yscale("log")
        ax.set_xlabel("N")
        ax.set_ylabel("Rel. amount of low-energy strength")
        # ax.set_ylim([1e-1, 4e0])
        ax.legend()
    
    plt.show()


if __name__ == "__main__":
    directory_base = "kshell_output.tmp"
    while True:
        directory = input("Choose kshell output directory: ")
        directory = f"{directory_base}/run_{directory}/"
        if os.path.isdir(directory):
            break
        else:
            print(f"'{directory}' is not a valid directory! Valid directories are:")
            print(os.listdir(directory_base))
    
    while True:
        filter_ = input("Choose element (atomic mass or name): ")
        
        if (filter_ == "") or (filter_ == "all"):
            filter_ = None
            break
        
        if any([filter_ in elem for elem in os.listdir(directory)]):
            try:
                filter_ = f"{ksutil.atomic_numbers[filter_]}_{filter_}"
            except KeyError:
                filter_ = f"{ksutil.atomic_numbers_reversed[int(filter_)]}_{filter_}"
            break
        else:
            print("Valid elements are:")
            print(os.listdir(directory))
    
    # q = Recreate(directory=directory)
    # q.recreate_figure_6(filter_=filter_)
    recreate(directory, filter_)
    # level_plot(directory=directory, filter_=filter_)

    # try:
    #     isotope_name = sys.argv[1]
    # except IndexError:
    #     isotope_name = "S24"
    # q.recreate_figure_6(filter=None)
    # q.plot_gsf(isotope_name)


    pass