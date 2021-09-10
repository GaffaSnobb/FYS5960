import sys, os
# from fractions import Fraction
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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

class Recreate:
    def __init__(self, directory):
        self.bin_width = 0.2
        self.E_max = 30
        self.Ex_min = 0 # Lower limit for emitted gamma energy [MeV].
        self.Ex_max = 30 # Upper limit for emitted gamma energy [MeV].
        n_bins = int(np.ceil(self.E_max/self.bin_width))
        E_max_adjusted = self.bin_width*n_bins
        bins = np.linspace(0, E_max_adjusted, n_bins + 1)
        self.bins_middle = (bins[0: -1] + bins[1:])/2

        self.all_fnames = {}

        # self.directory = "kshell_output/run_1/"
        self.directory = directory
        for element in sorted(os.listdir(self.directory)):
            """
            List all content in self.directory.
            """
            if os.path.isdir(self.directory + element):
                """
                If element is a directory, enter it to find data files.
                """
                self.all_fnames[element] = []    # Create blank entry in dict for current element.
                for isotope in os.listdir(self.directory + element):
                    """
                    List all content in the element directory.
                    """
                    if isotope.startswith("summary"):
                        """
                        Extract summary data files.
                        """
                        try:
                            """
                            Example: O16.
                            """
                            n_neutrons = int(isotope[9:11])
                        except ValueError:
                            """
                            Example: Ne20.
                            """
                            n_neutrons = int(isotope[10:12])

                        n_neutrons -= ksutil.atomic_numbers[element.split("_")[1]]
                        
                        self.all_fnames[element].append([element + "/" + isotope, n_neutrons])
        
        for key in self.all_fnames:
            """
            Sort each list in the dict by the number of neutrons.
            """
            self.all_fnames[key].sort(key=lambda tup: tup[1])   # Why not do this when directory is listed?

    def plot_gsf(self, isotope_name):
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
        res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with Jørgen.

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

    def recreate_figure_6(self, filter_=None):
        """
        Recreate the figure from Jørgens article.
        """
        fig, ax = plt.subplots()

        # for fnames in self.fnames_combined:
        for key in self.all_fnames:
            """
            Loop over all elements (grunnstoff).
            """
            fnames = self.all_fnames[key]   # For compatibility with old code.
            if filter_ is not None:
                if key.split("_")[1] not in filter_:
                    """
                    Skip elements not in filter_.
                    """
                    continue
            
            ratios = [] # Reset ratio for every new element.
            for i in range(len(fnames)):

                """
                Loop over all isotopes per element.
                """
                print(f"{fnames[i][0]=}")

                try:
                    res = ksutil.loadtxt(self.directory + fnames[i][0])[0]
                except FileNotFoundError:
                    print(f"File {fnames[i][0]} skipped! File not found.")
                    ratios.append(None) # Maintain correct list length for plotting.
                    continue

                Jpi_list = ksutil.create_jpi_list(
                    spins = res.levels[:, 1],
                    parities = res.levels[:, 2]
                )
                E_gs = res.levels[0, 0]

                try:
                    res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with Jørgen.
                except IndexError:
                    print(f"File {fnames[i][0]} skipped! Too few / no energy levels are present in this data file.")
                    ratios.append(None) # Maintain correct list length for plotting.
                    continue
                
                try:
                    gsf = ksutil.strength_function_average(
                        levels = res.levels,
                        transitions = res.transitions,
                        Jpi_list = Jpi_list,
                        bin_width = self.bin_width,
                        Ex_min = self.Ex_min,    # [MeV].
                        Ex_max = self.Ex_max,    # [MeV].
                        multipole_type = "M1"
                    )
                except IndexError:
                    print(f"File {fnames[i][0]} skipped! That unknown index out of bounds error in ksutil.")
                    ratios.append(None)
                    continue

                # Sum gsf for low and high energy range and take the ratio.
                bin_slice = self.bins_middle[0:len(gsf)]
                low_idx = (bin_slice <= 2)
                # high_idx = (bin_slice <= 6) == (2 <= bin_slice)
                high_idx = (bin_slice <= 100) == (2 <= bin_slice)
                low = np.sum(gsf[low_idx])
                high = np.sum(gsf[high_idx])
                low_high_ratio = low/high
                ratios.append(low_high_ratio)

            if all(elem is None for elem in ratios):
                """
                Skip current element if no ratios are calculated.
                """
                continue
            
            label = fnames[0][0][:fnames[0][0].index("/")]
            ax.plot([n_neutrons for _, n_neutrons in fnames], ratios, "--.", label=label)            
            ax.set_yscale("log")
            ax.set_xlabel("N")
            ax.set_ylabel("Rel. amount of low-energy strength")
            ax.set_ylim([1e-1, 4e0])
            ax.legend()
        
        plt.show()

def level_plot(directory, filter_):

    data_files = ksutil.loadtxt(path=directory, is_directory=True, filter_=filter_)
    print(data_files[6].fname_summary)
    print(data_files[6].levels)
    # line_width = 0.4
    # x_scope = range(len(data_files[1].E_x))
    # for i in range(len(data_files)):
    #     try:
    #         n = len(data_files[i].E_x)
    #     except TypeError:
    #         """
    #         E_x might be None.
    #         """
    #         continue
        
    #     for j in range(n):
    #         plt.hlines(
    #             y = data_files[i].E_x[j],
    #             xmin = x_scope[j] - line_width,
    #             xmax = x_scope[j] + line_width,
    #         )
    
    # plt.legend(loc="lower right")
    # plt.ylabel("MeV")
    # plt.show()

if __name__ == "__main__":
    directory_base = "kshell_output"
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
    
    q = Recreate(directory=directory)
    q.recreate_figure_6(filter_=filter_)
    # level_plot(directory=directory, filter_=filter_)

    # try:
    #     isotope_name = sys.argv[1]
    # except IndexError:
    #     isotope_name = "S24"
    # q.recreate_figure_6(filter=None)
    # q.plot_gsf(isotope_name)


    pass