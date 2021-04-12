import sys, os
from fractions import Fraction
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../../../../kshell_public/bin/')
import shellmodelutilities as smutil


atomic_numbers = {"oxygen": 8, "fluorine": 9, "neon": 10, "sodium": 11}


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


class Recreate:
    def __init__(self):
        self.bin_width = 0.2
        self.E_max = 30
        self.Ex_min = 0 # Lower limit for emitted gamma energy [MeV].
        self.Ex_max = 30 # Upper limit for emitted gamma energy [MeV].
        n_bins = int(np.ceil(self.E_max/self.bin_width))
        E_max_adjusted = self.bin_width*n_bins
        bins = np.linspace(0, E_max_adjusted, n_bins + 1)
        self.bins_middle = (bins[0: -1] + bins[1:])/2

        self.all_fnames = {}

        self.directory = "kshell_output/run_4/"
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

                        n_neutrons -= atomic_numbers[element.split("_")[1]]
                        
                        self.all_fnames[element].append([element + "/" + isotope, n_neutrons])
        
        for key in self.all_fnames:
            """
            Sort each list in the dict by the number of neutrons.
            """
            self.all_fnames[key].sort(key=lambda tup: tup[1])
        
        # self.fnames_oxygen_01 = [
        #     ["08_oxygen_01/summary_O16_usda.txt", 16], ["08_oxygen_01/summary_O17_usda.txt", 17],
        #     ["08_oxygen_01/summary_O18_usda.txt", 18], ["08_oxygen_01/summary_O19_usda.txt", 19],
        #     ["08_oxygen_01/summary_O20_usda.txt", 20], ["08_oxygen_01/summary_O21_usda.txt", 21],
        #     ["08_oxygen_01/summary_O22_usda.txt", 22], ["08_oxygen_01/summary_O23_usda.txt", 23],
        #     ["08_oxygen_01/summary_O24_usda.txt", 24], ["08_oxygen_01/summary_O25_usda.txt", 25],
        #     ["08_oxygen_01/summary_O26_usda.txt", 26],# ["08_oxygen_01/summary_O27_usda.txt", 27],
        #     #["08_oxygen_01/summary_O28_usda.txt", 28]
        # ]
        # self.fnames_oxygen_02 = [
        #     ["08_oxygen_02/summary_O16_usda.txt", 16], ["08_oxygen_02/summary_O17_usda.txt", 17],
        #     ["08_oxygen_02/summary_O18_usda.txt", 18], ["08_oxygen_02/summary_O19_usda.txt", 19],
        #     ["08_oxygen_02/summary_O20_usda.txt", 20], ["08_oxygen_02/summary_O21_usda.txt", 21],
        #     ["08_oxygen_02/summary_O22_usda.txt", 22], ["08_oxygen_02/summary_O23_usda.txt", 23],
        #     ["08_oxygen_02/summary_O24_usda.txt", 24], ["08_oxygen_02/summary_O25_usda.txt", 25],
        #     ["08_oxygen_02/summary_O26_usda.txt", 26],# ["08_oxygen_02/summary_O27_usda.txt", 27],
        #     #["08_oxygen_02/summary_O28_usda.txt", 28]
        # ]
        # self.fnames_oxygen_03 = [
        #     ["08_oxygen_03/summary_O16_usda.txt", 16], ["08_oxygen_03/summary_O17_usda.txt", 17],
        #     ["08_oxygen_03/summary_O18_usda.txt", 18], ["08_oxygen_03/summary_O19_usda.txt", 19],
        #     ["08_oxygen_03/summary_O20_usda.txt", 20], ["08_oxygen_03/summary_O21_usda.txt", 21],
        #     ["08_oxygen_03/summary_O22_usda.txt", 22], ["08_oxygen_03/summary_O23_usda.txt", 23],
        #     ["08_oxygen_03/summary_O24_usda.txt", 24], ["08_oxygen_03/summary_O25_usda.txt", 25],
        #     ["08_oxygen_03/summary_O26_usda.txt", 26],# ["08_oxygen_03/summary_O27_usda.txt", 27],
        #     #["08_oxygen_03/summary_O28_usda.txt", 28]
        # ]
        # self.fnames_oxygen = [
        #     ["08_oxygen/summary_O16_usda.txt", 16], ["08_oxygen/summary_O17_usda.txt", 17],
        #     ["08_oxygen/summary_O18_usda.txt", 18], ["08_oxygen/summary_O19_usda.txt", 19],
        #     ["08_oxygen/summary_O20_usda.txt", 20], ["08_oxygen/summary_O21_usda.txt", 21],
        #     ["08_oxygen/summary_O22_usda.txt", 22], ["08_oxygen/summary_O23_usda.txt", 23],
        #     ["08_oxygen/summary_O24_usda.txt", 24], ["08_oxygen/summary_O25_usda.txt", 25],
        #     ["08_oxygen/summary_O26_usda.txt", 26],# ["08_oxygen/summary_O27_usda.txt", 27],
        #     #["08_oxygen/summary_O28_usda.txt", 28]
        # ]
        # self.fnames_fluorine = [
        #     ["09_fluorine/summary_F17_usda.txt", 17], ["09_fluorine/summary_F18_usda.txt", 18],
        #     ["09_fluorine/summary_F19_usda.txt", 19], ["09_fluorine/summary_F20_usda.txt", 20],
        #     ["09_fluorine/summary_F21_usda.txt", 21], ["09_fluorine/summary_F22_usda.txt", 22],
        #     ["09_fluorine/summary_F23_usda.txt", 23], ["09_fluorine/summary_F24_usda.txt", 24],
        #     ["09_fluorine/summary_F25_usda.txt", 25], ["09_fluorine/summary_F26_usda.txt", 26],
        #     ["09_fluorine/summary_F27_usda.txt", 27], ["09_fluorine/summary_F28_usda.txt", 28],
        #     ["09_fluorine/summary_F29_usda.txt", 29]
        # ]
        # self.fnames_neon = [
        #     ["10_neon/summary_Ne18_usda.txt", 18], ["10_neon/summary_Ne19_usda.txt", 19],
        #     ["10_neon/summary_Ne20_usda.txt", 20], ["10_neon/summary_Ne21_usda.txt", 21],
        #     ["10_neon/summary_Ne22_usda.txt", 22], ["10_neon/summary_Ne23_usda.txt", 23],
        #     ["10_neon/summary_Ne24_usda.txt", 24], ["10_neon/summary_Ne25_usda.txt", 25],
        #     ["10_neon/summary_Ne26_usda.txt", 26], ["10_neon/summary_Ne27_usda.txt", 27],
        #     ["10_neon/summary_Ne28_usda.txt", 28], ["10_neon/summary_Ne29_usda.txt", 29],
        #     ["10_neon/summary_Ne30_usda.txt", 30]
        # ]
        # self.fnames_neon_01 = [
        #     ["10_neon_01/summary_Ne18_usda.txt", 18], ["10_neon_01/summary_Ne19_usda.txt", 19],
        #     ["10_neon_01/summary_Ne20_usda.txt", 20], ["10_neon_01/summary_Ne21_usda.txt", 21],
        #     ["10_neon_01/summary_Ne22_usda.txt", 22], ["10_neon_01/summary_Ne23_usda.txt", 23],
        #     ["10_neon_01/summary_Ne24_usda.txt", 24], ["10_neon_01/summary_Ne25_usda.txt", 25],
        #     ["10_neon_01/summary_Ne26_usda.txt", 26], ["10_neon_01/summary_Ne27_usda.txt", 27],
        #     ["10_neon_01/summary_Ne28_usda.txt", 28], ["10_neon_01/summary_Ne29_usda.txt", 29],
        #     ["10_neon_01/summary_Ne30_usda.txt", 30]
        # ]
        # self.fnames_neon_02 = [
        #     ["10_neon_02/summary_Ne18_usda.txt", 18], ["10_neon_02/summary_Ne19_usda.txt", 19],
        #     ["10_neon_02/summary_Ne20_usda.txt", 20], ["10_neon_02/summary_Ne21_usda.txt", 21],
        #     ["10_neon_02/summary_Ne22_usda.txt", 22], ["10_neon_02/summary_Ne23_usda.txt", 23],
        #     ["10_neon_02/summary_Ne24_usda.txt", 24], ["10_neon_02/summary_Ne25_usda.txt", 25],
        #     ["10_neon_02/summary_Ne26_usda.txt", 26], ["10_neon_02/summary_Ne27_usda.txt", 27],
        #     ["10_neon_02/summary_Ne28_usda.txt", 28], ["10_neon_02/summary_Ne29_usda.txt", 29],
        #     ["10_neon_02/summary_Ne30_usda.txt", 30]
        # ]
        # self.fnames_neon_03 = [
        #     ["10_neon_03/summary_Ne18_usda.txt", 18], ["10_neon_03/summary_Ne19_usda.txt", 19],
        #     ["10_neon_03/summary_Ne20_usda.txt", 20], ["10_neon_03/summary_Ne21_usda.txt", 21],
        #     ["10_neon_03/summary_Ne22_usda.txt", 22], ["10_neon_03/summary_Ne23_usda.txt", 23],
        #     ["10_neon_03/summary_Ne24_usda.txt", 24], ["10_neon_03/summary_Ne25_usda.txt", 25],
        #     ["10_neon_03/summary_Ne26_usda.txt", 26], ["10_neon_03/summary_Ne27_usda.txt", 27],
        #     ["10_neon_03/summary_Ne28_usda.txt", 28], ["10_neon_03/summary_Ne29_usda.txt", 29],
        #     ["10_neon_03/summary_Ne30_usda.txt", 30]
        # ]
        # self.fnames_sodium = [
        #     ["11_sodium/summary_Na19_usda.txt", 19], ["11_sodium/summary_Na20_usda.txt", 20],
        #     ["11_sodium/summary_Na21_usda.txt", 21], ["11_sodium/summary_Na22_usda.txt", 22],
        #     ["11_sodium/summary_Na23_usda.txt", 23], ["11_sodium/summary_Na24_usda.txt", 24],
        #     ["11_sodium/summary_Na25_usda.txt", 25], ["11_sodium/summary_Na26_usda.txt", 26],
        #     ["11_sodium/summary_Na27_usda.txt", 27], ["11_sodium/summary_Na28_usda.txt", 28],
        #     ["11_sodium/summary_Na29_usda.txt", 29], ["11_sodium/summary_Na30_usda.txt", 30],
        #     ["11_sodium/summary_Na31_usda.txt", 31]
        # ]
        # self.fnames_magnesium = [
        #     ["12_magnesium/summary_Mg20_usda.txt", 20], ["12_magnesium/summary_Mg21_usda.txt", 21],
        #     ["12_magnesium/summary_Mg22_usda.txt", 22], ["12_magnesium/summary_Mg23_usda.txt", 23],
        #     ["12_magnesium/summary_Mg24_usda.txt", 24], ["12_magnesium/summary_Mg25_usda.txt", 25],
        #     ["12_magnesium/summary_Mg26_usda.txt", 26], ["12_magnesium/summary_Mg27_usda.txt", 27],
        #     ["12_magnesium/summary_Mg28_usda.txt", 28], ["12_magnesium/summary_Mg29_usda.txt", 29],
        #     ["12_magnesium/summary_Mg30_usda.txt", 30], ["12_magnesium/summary_Mg31_usda.txt", 31],
        #     # ["12_magnesium/summary_Mg32_usda.txt", 32]
        # ]
        # self.fnames_aluminium = [
        #     ["13_aluminium/summary_Al21_usda.txt", 21], ["13_aluminium/summary_Al22_usda.txt", 22],
        #     ["13_aluminium/summary_Al23_usda.txt", 23], ["13_aluminium/summary_Al24_usda.txt", 24],
        #     ["13_aluminium/summary_Al25_usda.txt", 25], ["13_aluminium/summary_Al26_usda.txt", 26],
        #     ["13_aluminium/summary_Al27_usda.txt", 27], ["13_aluminium/summary_Al28_usda.txt", 28],
        #     ["13_aluminium/summary_Al29_usda.txt", 29], ["13_aluminium/summary_Al30_usda.txt", 30],
        #     ["13_aluminium/summary_Al31_usda.txt", 31], ["13_aluminium/summary_Al32_usda.txt", 32],
        #     ["13_aluminium/summary_Al33_usda.txt", 33]
        # ]
        # self.fnames_silicon = [
        #     ["14_silicon/summary_Si22_usda.txt", 22], ["14_silicon/summary_Si23_usda.txt", 23],
        #     ["14_silicon/summary_Si24_usda.txt", 24], ["14_silicon/summary_Si25_usda.txt", 25],
        #     ["14_silicon/summary_Si26_usda.txt", 26], ["14_silicon/summary_Si27_usda.txt", 27],
        #     ["14_silicon/summary_Si28_usda.txt", 28], ["14_silicon/summary_Si29_usda.txt", 29],
        #     ["14_silicon/summary_Si30_usda.txt", 30], ["14_silicon/summary_Si31_usda.txt", 31],
        #     ["14_silicon/summary_Si32_usda.txt", 32], ["14_silicon/summary_Si33_usda.txt", 33],
        #     ["14_silicon/summary_Si34_usda.txt", 34]
        # ]
        # self.fnames_phosphorus = [
        #     ["15_phosphorus/summary_P23_usda.txt", 23], ["15_phosphorus/summary_P24_usda.txt", 24],
        #     ["15_phosphorus/summary_P25_usda.txt", 25], ["15_phosphorus/summary_P26_usda.txt", 26],
        #     ["15_phosphorus/summary_P27_usda.txt", 27], ["15_phosphorus/summary_P28_usda.txt", 28],
        #     ["15_phosphorus/summary_P29_usda.txt", 29], ["15_phosphorus/summary_P30_usda.txt", 30],
        #     ["15_phosphorus/summary_P31_usda.txt", 31], ["15_phosphorus/summary_P32_usda.txt", 32],
        #     ["15_phosphorus/summary_P33_usda.txt", 33], ["15_phosphorus/summary_P34_usda.txt", 34],
        #     ["15_phosphorus/summary_P35_usda.txt", 35]
        # ]
        # self.fnames_sulfur = [
        #     ["16_sulfur/summary_S24_usda.txt", 24], ["16_sulfur/summary_S25_usda.txt", 25],
        #     ["16_sulfur/summary_S26_usda.txt", 26], ["16_sulfur/summary_S27_usda.txt", 27],
        #     ["16_sulfur/summary_S28_usda.txt", 28], ["16_sulfur/summary_S29_usda.txt", 29],
        #     ["16_sulfur/summary_S30_usda.txt", 30], ["16_sulfur/summary_S31_usda.txt", 31],
        #     ["16_sulfur/summary_S32_usda.txt", 32], ["16_sulfur/summary_S33_usda.txt", 33],
        #     ["16_sulfur/summary_S34_usda.txt", 34], ["16_sulfur/summary_S35_usda.txt", 35],
        #     ["16_sulfur/summary_S36_usda.txt", 36]
        # ]
        # self.fnames_chlorine = [
        #     ["17_chlorine/summary_Cl26_usda.txt", 26], ["17_chlorine/summary_Cl27_usda.txt", 27],
        #     ["17_chlorine/summary_Cl28_usda.txt", 28], ["17_chlorine/summary_Cl29_usda.txt", 29],
        #     ["17_chlorine/summary_Cl30_usda.txt", 30], ["17_chlorine/summary_Cl31_usda.txt", 31],
        #     ["17_chlorine/summary_Cl32_usda.txt", 32], ["17_chlorine/summary_Cl33_usda.txt", 33],
        #     ["17_chlorine/summary_Cl34_usda.txt", 34], ["17_chlorine/summary_Cl35_usda.txt", 35],
        #     ["17_chlorine/summary_Cl36_usda.txt", 36], ["17_chlorine/summary_Cl37_usda.txt", 37]
        # ]
        # self.fnames_argon = [
        #     ["18_argon/summary_Ar26_usda.txt", 26], ["18_argon/summary_Ar27_usda.txt", 27],
        #     ["18_argon/summary_Ar28_usda.txt", 28], ["18_argon/summary_Ar29_usda.txt", 29],
        #     ["18_argon/summary_Ar30_usda.txt", 30], ["18_argon/summary_Ar31_usda.txt", 31],
        #     ["18_argon/summary_Ar32_usda.txt", 32], ["18_argon/summary_Ar33_usda.txt", 33],
        #     ["18_argon/summary_Ar34_usda.txt", 34], ["18_argon/summary_Ar35_usda.txt", 35],
        #     ["18_argon/summary_Ar36_usda.txt", 36], ["18_argon/summary_Ar37_usda.txt", 37],
        #     ["18_argon/summary_Ar38_usda.txt", 38],
        # ]
        # self.fnames_oxygen_01 = [[fname, mass_number - 8 ] for fname, mass_number in self.fnames_oxygen_01]
        # self.fnames_oxygen_02 = [[fname, mass_number - 8 ] for fname, mass_number in self.fnames_oxygen_02]
        # self.fnames_oxygen_03 = [[fname, mass_number - 8 ] for fname, mass_number in self.fnames_oxygen_03]
        # self.fnames_neon_01   = [[fname, mass_number - 10] for fname, mass_number in self.fnames_neon_01]
        # self.fnames_neon_02   = [[fname, mass_number - 10] for fname, mass_number in self.fnames_neon_02]
        # self.fnames_neon_03   = [[fname, mass_number - 10] for fname, mass_number in self.fnames_neon_03]
        
        # self.fnames_oxygen     = [[fname, mass_number - 8 ] for fname, mass_number in self.fnames_oxygen]
        # self.fnames_fluorine   = [[fname, mass_number - 9 ] for fname, mass_number in self.fnames_fluorine]
        # self.fnames_neon       = [[fname, mass_number - 10] for fname, mass_number in self.fnames_neon]
        # self.fnames_sodium     = [[fname, mass_number - 11] for fname, mass_number in self.fnames_sodium]
        # self.fnames_magnesium  = [[fname, mass_number - 12] for fname, mass_number in self.fnames_magnesium]
        # self.fnames_aluminium  = [[fname, mass_number - 13] for fname, mass_number in self.fnames_aluminium]
        # self.fnames_silicon    = [[fname, mass_number - 14] for fname, mass_number in self.fnames_silicon]
        # self.fnames_phosphorus = [[fname, mass_number - 15] for fname, mass_number in self.fnames_phosphorus]
        # self.fnames_sulfur     = [[fname, mass_number - 16] for fname, mass_number in self.fnames_sulfur]
        # self.fnames_chlorine   = [[fname, mass_number - 17] for fname, mass_number in self.fnames_chlorine]
        # self.fnames_argon      = [[fname, mass_number - 18] for fname, mass_number in self.fnames_argon]
        
        # self.fnames_combined = [
        #     self.fnames_oxygen_01, self.fnames_oxygen_02, self.fnames_oxygen_03,
        #     self.fnames_neon_01, self.fnames_neon_02, self.fnames_neon_03
        # ]
        # self.fnames_combined = [
        #     self.fnames_oxygen,   self.fnames_fluorine,   self.fnames_neon,
        #     self.fnames_sodium,   self.fnames_magnesium,  self.fnames_aluminium,
        #     self.fnames_silicon,  self.fnames_phosphorus, self.fnames_sulfur,
        #     self.fnames_chlorine, self.fnames_argon
        # ]

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

        res = ReadKshellOutput()
        res.read_kshell_output(fname = self.directory + fname)

        fig, ax = plt.subplots()

        Jpi_list = create_jpi_list(res.levels[:, 1], None)
        E_gs = res.levels[0, 0]
        res.transitions[:, 2] += E_gs   # Add ground state energy for compatibility with Jørgen.

        gsf = smutil.strength_function_average(
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


    def recreate_figure_6(self, filter=None):
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
            if filter is not None:
                if key.split("_")[1] not in filter:
                    continue
            ratios = [] # Reset ratio for every new element.
            for i in range(len(fnames)):

                """
                Loop over all isotopes per element.
                """
                print(f"{fnames[i][0]=}")
                res = ReadKshellOutput()

                try:
                    res.read_kshell_output(fname = self.directory + fnames[i][0])
                except FileNotFoundError:
                    print(f"File {fnames[i][0]} skipped! File not found.")
                    ratios.append(None) # Maintain correct list length for plotting.
                    continue

                Jpi_list = create_jpi_list(
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

                gsf = smutil.strength_function_average(
                    levels = res.levels,
                    transitions = res.transitions,
                    Jpi_list = Jpi_list,
                    bin_width = self.bin_width,
                    Ex_min = self.Ex_min,    # [MeV].
                    Ex_max = self.Ex_max,    # [MeV].
                    multipole_type = "M1"
                )

                # Sum gsf for low and high energy range and take the ratio.
                bin_slice = self.bins_middle[0:len(gsf)]
                low_idx = (bin_slice <= 2)
                high_idx = (bin_slice <= 6) == (2 <= bin_slice)
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
            ax.legend()
        
        plt.show()


if __name__ == "__main__":
    test_file_read()
    
    try:
        isotope_name = sys.argv[1]
    except IndexError:
        isotope_name = "S24"
    
    q = Recreate()
    q.recreate_figure_6()
    # q.plot_gsf(isotope_name)


    pass