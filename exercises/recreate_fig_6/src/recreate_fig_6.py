import sys
import numpy as np
import matplotlib.pyplot as plt

def level_density_lol(bin_size):
    """
    Read level data from file.  Calculate the level density for input
    bin size.

    Parameters
    ----------
    bin_size : float, int
        The bin size.  Must be same unit as level data from file.

    Returns
    -------
    counts : numpy.ndarray
        The number of levels within +- bin_size for all levels.

    counts/bin_size : numpy.ndarray
        The level density.

    levels : numpy.ndarray
        The raw level data.
    """
    fname = "levels_56Fe.txt"
    levels = np.loadtxt(fname)
    n_levels = len(levels)
    counts = np.zeros(n_levels)

    for i in range(n_levels):
        for j in range(n_levels):
            if (levels[j] >= levels[i] - bin_size) and (levels[j] <= levels[i] + bin_size):
                counts[i] += 1

    return counts, counts/bin_size, levels


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
    counts/bin_sizes : numpy.ndarray
        The level density.
    """
    N = len(levels)
    counts = []
    bin_levels = []
    current_bin = 0  # Initial.
    bins = np.arange(0, levels[-1] + bin_size, bin_size)
    counts = np.zeros(len(bins))
    # for i in range(N):
    #     for j in range(N):
    #         if (levels[j] >= levels[i] - bin_size/2) and (levels[j] <= levels[i] + bin_size/2):
    #             counts[i] += 1

    # while current_bin < levels[-1]:
    #     bin_levels.append(current_bin)
    #     counts.append(np.sum(  current_bin <= levels[levels < current_bin + bin_size]  ))
    #     current_bin += bin_size

    for i in range(len(bins) - 1):
        counts[i] = np.sum(bins[i] <= levels[levels < bins[i + 1]])

    return (counts/bin_size)[:-1], bins[1:]


    # counts = np.array(counts)
    # return counts/bin_size, np.array(bin_levels)




def read_kshell_output(fname):
    """
    Read energy level data, transition probabilities and transition
    strengths from KSHELL output files.

    Parameters
    ----------
    fname : string
        Filename of KSHELL output file.

    Returns
    -------
    E_x : list
        List of energy levels.

    B_M1 : list
        Nested list, [E_x, b_m1, E_gamma].
    """

    E_x = []    
    B_M1 = []

    def load_energy_levels(infile):
        for _ in range(3): infile.readline()
        for line in infile:
            try:
                tmp = line.split()
                E_x.append(float(tmp[6]))
            except IndexError:
                """
                End of energies.
                """
                break

    def load_m1_probabilities(infile):
        for _ in range(2): infile.readline()
        for line in infile:
            try:
                tmp = line.split()
                if (tmp[1][-1] != ")") and (tmp[3][-1] != ")"):
                    """
                    Example: 2+(11) 18.393 2+(10) 17.791 0.602 0.1(  0.0) 0.1( 0.0)
                    """
                    E_gamma = abs(float(tmp[3]) - float(tmp[1]))
                    E_i = float(tmp[1])
                    b_m1 = float(tmp[5][:-1])
                    B_M1.append([E_i, b_m1, E_gamma])

                elif (tmp[1][-1] != ")") and (tmp[3][-1] == ")"):
                    """
                    Example: 2+(10) 17.791 2+( 1) 5.172 12.619 0.006( 0.00) 0.006( 0.00)
                    """
                    E_gamma = abs(float(tmp[4]) - float(tmp[1]))
                    E_i = float(tmp[1])
                    b_m1 = float(tmp[6][:-1])
                    B_M1.append([E_i, b_m1, E_gamma])
                
                elif (tmp[1][-1] == ")") and (tmp[4][-1] != ")"):
                    """
                    Example: 3+( 8) 19.503 2+(11) 18.393 1.111 0.000( 0.00) 0.000( 0.00)
                    """
                    E_gamma = abs(float(tmp[4]) - float(tmp[2]))
                    E_i = float(tmp[2])
                    b_m1 = float(tmp[6][:-1])
                    B_M1.append([E_i, b_m1, E_gamma])

                elif (tmp[1][-1] == ")") and (tmp[4][-1] == ")"):
                    """
                    Example: 1+( 7) 19.408 2+( 9) 16.111 3.297 0.005( 0.00) 0.003( 0.00)
                    """
                    E_gamma = abs(float(tmp[5]) - float(tmp[2]))
                    E_i = float(tmp[2])
                    b_m1 = float(tmp[7][:-1])
                    B_M1.append([E_i, b_m1, E_gamma])
                else:
                    print("WARNING: Structure not accounted for! Exiting...")
                    print(line)
                    sys.exit()

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

    return E_x, B_M1


def test_file_read():
    """
    Test that 'read_kshell_output' successfully reads output from kshell.

    Raises
    ------
    AssertionError:
        If the read values are not exactly equal to the expected values.
    """
    fname = "test_text_file.txt"
    E_x, B_M1 = read_kshell_output(fname)
    E_x_expected = [0.0, 8.016]
    B_M1_expected = [
        [5.172, 20.5, 5.172],
        [17.791, 0.0, 17.791],
        [19.408, 5.7, 19.408 - 17.791],
        [18.393, 0.1, 18.393 - 17.791]
    ]
    
    for calculated, expected in zip(E_x, E_x_expected):
        msg = f"Error in E_x. Expected: {expected}, got: {calculated}."
        assert calculated == expected, msg

    for calculated, expected in zip(B_M1, B_M1_expected):
        msg = f"Error in B_M1. Expected: {expected}, got: {calculated}."
        success = (calculated[0] == expected[0]) and (calculated[1] == expected[1]) and (calculated[2] == expected[2])
        assert success, msg


if __name__ == "__main__":
    test_file_read()
    # levels_1, B_M1 = read_kshell_output("kshell_output/summary_S24_usda.txt")
    bin_size_1 = 0.2    # [MeV].
    bin_size_2 = 0.5    # [MeV].

    levels_2 = np.loadtxt("levels_56Fe.txt")
    densities_1, bin_levels_1 = level_density(levels_2, bin_size_1*1000)
    densities_2, bin_levels_2 = level_density(levels_2, bin_size_2*1000)

    print(bin_levels_2)
    
    plt.step(bin_levels_1, densities_1, label=f"{bin_size_1=}")
    plt.step(bin_levels_2, densities_2, label=f"{bin_size_2=}")
    plt.yscale("log")
    plt.legend()
    plt.savefig("level_density.png", dpi=300)
    plt.show()



    pass