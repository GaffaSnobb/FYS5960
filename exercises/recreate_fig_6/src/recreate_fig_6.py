import sys
from fractions import Fraction
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../../../../kshell_public/bin/')
import shellmodelutilities as smutil


def level_density(levels, bin_size):
    """
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

class ReadKshellOutput:
    """
    Implemented as class just to avoid returning a bunch of values.
    Access instance attributes instead.
    """
    def read_kshell_output(self, fname):
        """
        Read energy level data, transition probabilities and transition
        strengths from KSHELL output files.

        Parameters
        ----------
        fname : string
            Filename of KSHELL output file.

        Returns
        -------
        E_x : numpy.ndarray
            1D array of energy levels.

        B_M1 : numpy.ndarray
            Nx3 matrix with rows [E_x, b_m1, E_gamma].

        Raises
        ------
        RuntimeError
            If the KSHELL file has unexpected structure / syntax.
        """

        self.E_x = []    
        self.B_M1 = []
        self.levels = [] # [Ei, 2*Ji, parity].
        self.transitions = []   # [2J_f, p_i, E_f, 2J_i, p_i, E_i, E_gamma, B(.., i->f)].

        def load_energy_levels(infile):
            for _ in range(3): infile.readline()
            for line in infile:
                try:
                    tmp = line.split()
                    self.E_x.append(float(tmp[6]))
                    parity = 1 if tmp[2] == "+" else -1
                    self.levels.append([float(tmp[5]), 2*float(Fraction(tmp[1])), parity])
                except IndexError:
                    """
                    End of energies.
                    """
                    break

        def load_m1_probabilities(infile):
            """
            """
            for _ in range(2): infile.readline()
            for line in infile:
                try:
                    """
                    Example of possible lines in file:
                    J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
                    2+(11) 18.393 2+(10) 17.791 0.602 0.1(  0.0) 0.1( 0.0)
                    3/2+( 1) 0.072 5/2+( 1) 0.000 0.071 0.127( 0.07) 0.084( 0.05)
                    2+(10) 17.791 2+( 1) 5.172 12.619 0.006( 0.00) 0.006( 0.00)
                    3+( 8) 19.503 2+(11) 18.393 1.111 0.000( 0.00) 0.000( 0.00)
                    1+( 7) 19.408 2+( 9) 16.111 3.297 0.005( 0.00) 0.003( 0.00)
                    """
                    tmp = line.split()
                    
                    # Location of initial parity is common for all cases.
                    parity_idx = tmp[0].index("(") - 1 # Find index of initial parity.
                    p_i = 1 if tmp[0][parity_idx] == "+" else -1
                    
                    # Location of initial spin is common for all cases.
                    J_i = float(Fraction(tmp[0][:parity_idx]))
                    
                    if (tmp[1][-1] != ")") and (tmp[3][-1] != ")"):
                        """
                        Example:
                        J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
                        2+(11) 18.393 2+(10) 17.791 0.602 0.1(  0.0) 0.1( 0.0)
                        """
                        E_gamma = float(tmp[4])
                        E_i = float(tmp[1])
                        b_m1 = float(tmp[5][:-1])
                        J_f = float(Fraction(tmp[2][:-5]))
                        E_f = float(tmp[3])
                        self.B_M1.append([E_i, b_m1, E_gamma])
                        self.transitions.append([2*J_f, p_i, E_f, 2*J_i, p_i, E_i, E_gamma, b_m1])

                    elif (tmp[1][-1] != ")") and (tmp[3][-1] == ")"):
                        """
                        Example:
                        J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
                        2+(10) 17.791 2+( 1) 5.172 12.619 0.006( 0.00) 0.006( 0.00)
                        """
                        E_gamma = float(tmp[5])
                        E_i = float(tmp[1])
                        b_m1 = float(tmp[6][:-1])
                        J_f = float(Fraction(tmp[2][:-2]))
                        E_f = float(tmp[4])
                        self.B_M1.append([E_i, b_m1, E_gamma])
                        self.transitions.append([2*J_f, p_i, E_f, 2*J_i, p_i, E_i, E_gamma, b_m1])
                    
                    elif (tmp[1][-1] == ")") and (tmp[4][-1] != ")"):
                        """
                        Example:
                        J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
                        3+( 8) 19.503 2+(11) 18.393 1.111 0.000( 0.00) 0.000( 0.00)
                        """
                        E_gamma = float(tmp[5])
                        E_i = float(tmp[2])
                        b_m1 = float(tmp[6][:-1])
                        J_f = float(Fraction(tmp[3][:-5]))
                        E_f = float(tmp[4])
                        self.B_M1.append([E_i, b_m1, E_gamma])
                        self.transitions.append([2*J_f, p_i, E_f, 2*J_i, p_i, E_i, E_gamma, b_m1])

                    elif (tmp[1][-1] == ")") and (tmp[4][-1] == ")"):
                        """
                        Example:
                        J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
                        1+( 7) 19.408 2+( 9) 16.111 3.297 0.005( 0.00) 0.003( 0.00)
                        """
                        E_gamma = float(tmp[6])
                        E_i = float(tmp[2])
                        b_m1 = float(tmp[7][:-1])
                        J_f = float(Fraction(tmp[3][:-2]))
                        E_f = float(tmp[5])
                        self.B_M1.append([E_i, b_m1, E_gamma])
                        self.transitions.append([2*J_f, p_i, E_f, 2*J_i, p_i, E_i, E_gamma, b_m1])
                    else:
                        msg = "WARNING: Structure not accounted for!"
                        msg += f"\n{line=}"
                        raise RuntimeError(msg)

                except IndexError:
                    """
                    End of probabilities.
                    """
                    break

        with open(fname, "r") as infile:
            for line in infile:
                tmp = line.split()
                try:
                    if tmp[0] == "Energy":
                        load_energy_levels(infile)
                    
                    elif tmp[0] == "B(E2)":
                        continue
                    
                    elif tmp[0] == "B(M1)":
                        load_m1_probabilities(infile)
                
                except IndexError:
                    """
                    Skip blank lines.
                    """
                    continue


        self.levels = np.array(self.levels)
        self.transitions = np.array(self.transitions)
        return np.array(self.E_x), np.array(self.B_M1)


def test_file_read():
    """
    Test that 'read_kshell_output' successfully reads output from kshell.

    Raises
    ------
    AssertionError
        If the read values are not exactly equal to the expected values.
    """
    res = ReadKshellOutput()
    E_x, B_M1 = res.read_kshell_output("test_text_file.txt")
    E_x_expected = [0.0, 8.016]
    B_M1_expected = [
        [5.172, 20.5, 5.172],
        [17.791, 0.0, 17.791],
        [19.408, 5.7, 1.617],
        [18.393, 0.1, 0.602]
    ]
    
    for calculated, expected in zip(E_x, E_x_expected):
        msg = f"Error in E_x. Expected: {expected}, got: {calculated}."
        assert calculated == expected, msg

    for calculated, expected in zip(B_M1, B_M1_expected):
        msg = f"Error in B_M1. Expected: {expected}, got: {calculated}."
        success = (calculated[0] == expected[0]) and (calculated[1] == expected[1]) and (calculated[2] == expected[2])
        assert success, msg


def mean_transition_strength(B_M1):
    """
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


def create_jpi_list(spins, parities):
    """
    Example list:
    [[1, +1], [3, +1], [5, +1], [7, +1], [9, +1], [11, +1], [13, +1]].
    Currently hard-coded for only positive parity states.

    Parameters
    ----------
    spins : numpy.ndarray
        Array of spins for each energy level.

    Returns
    -------
    spins_new : list
        A nested list of spins and parities [[spin, parity], ...] sorted
        with respect to the spin.
    """
    spins[spins < 0] = 0    # Discard negative entries.
    spins_new = []
    for elem in spins:
        if elem not in spins_new:
            """
            Extract each distinct value only once.
            """
            spins_new.append([elem])

    for i in range(len(spins_new)):
        """
        Add parity for each spin. [spin, parity].
        """
        spins_new[i].append(+1)
    
    return sorted(spins_new, key=lambda tup: tup[0])

if __name__ == "__main__":
    test_file_read()
    bin_width = 0.2
    E_max = 30
    Ex_min = 0 # Lower limit for emitted gamma energy [MeV].
    Ex_max = 20 # Upper limit for emitted gamma energy [MeV].
    n_bins = int(np.ceil(E_max/bin_width))
    E_max_adjusted = bin_width*n_bins
    bins = np.linspace(0, E_max_adjusted, n_bins + 1)
    bins_middle = (bins[0: -1] + bins[1:])/2

    fig, ax = plt.subplots()

    directory = "kshell_output/"
    fnames = [
        "summary_S26_usda.txt", "summary_S27_usda.txt", "summary_S28_usda.txt",
        "summary_S29_usda.txt", "summary_S30_usda.txt", "summary_S31_usda.txt",
        "summary_S32_usda.txt", "summary_S33_usda.txt", "summary_S34_usda.txt",
        "summary_S35_usda.txt",
    ]
    mass_numbers = np.array([26, 27, 28, 29, 30, 31, 32, 33, 34, 35])
    mass_numbers -= 16
    ratios = []

    for i in range(len(fnames)):
        res = ReadKshellOutput()
        _, _ = res.read_kshell_output(fname = directory + fnames[i])

        Jpi_list = create_jpi_list(
            spins = res.levels[:, 1],
            parities = res.levels[:, 2]
        )

        E_gs = res.levels[0, 0]
        res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with Jørgen.

        gsf = smutil.strength_function_average(
            levels = res.levels,
            transitions = res.transitions,
            Jpi_list = Jpi_list,
            bin_width = bin_width,
            Ex_min = Ex_min,    # [MeV].
            Ex_max = Ex_max,    # [MeV].
            multipole_type = "M1"
        )
        bin_slice = bins_middle[0:len(gsf)]
        low_idx = (bin_slice <= 2)
        high_idx = (bin_slice <= 6) == (2 <= bin_slice)
        low = np.sum(gsf[low_idx])
        high = np.sum(gsf[high_idx])

        low_high_ratio = low/high

        print(f"{low_high_ratio=}")
        ratios.append(low_high_ratio)

        # ax.plot(bin_slice, gsf, label=fnames[i][8:11])
        # ax.legend()
        # ax.set_xlabel(r"$E_{\gamma}$ [MeV]")
        # ax.set_ylabel(r"gsf [MeV$^{-3}$]")
    

    ax.plot(mass_numbers, ratios, "--.", label="S")
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel("Rel. amount of low-energy strength")
    ax.legend()
    plt.show()

    # print(f"{res.B_M1=}")
    # bin_size_1 = 0.2    # [MeV].
    # bin_size_2 = 0.5    # [MeV].

    # densities_1, bin_levels_1 = level_density(levels_1, bin_size_1)
    # densities_2, bin_levels_2 = level_density(levels_1, bin_size_2)

    # B_M1_mean, B_M1_mean_levels = mean_transition_strength(B_M1)

    # # plt.plot(B_M1_mean_levels, B_M1_mean, ".")
    # # plt.xlabel("Enegy level [MeV]")
    # # plt.ylabel(r"$\langle B(M1) \rangle$")
    # # # plt.yscale("log")
    # # plt.show()
    
    # plt.step(bin_levels_1, densities_1, label=f"{bin_size_1=}")
    # plt.step(bin_levels_2, densities_2, label=f"{bin_size_2=}")
    # plt.yscale("log")
    # plt.legend()
    # plt.savefig("level_density.png", dpi=300)
    # plt.show()



    pass