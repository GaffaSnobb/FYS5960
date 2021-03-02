import numpy as np

def level_density(bin_size):
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
    fname = ""
    levels = np.loadtxt(fname)
    n_levels = len(levels)
    counts = np.zeros(n_levels)

    for i in range(n_levels):
        for j in range(n_levels):
            if (levels[j] >= levels[i] - bin_size) and (levels[j] <= levels[i] + bin_size):
                counts[i] += 1

    return counts, counts/bin_size, levels


def read_kshell_output(fname="summary_S24_usda.txt"):

    E_x = []    
    B_M1 = []

    def load_energy_levels(infile):
        for _ in range(3): infile.readline()
        for line in infile:
            try:
                tmp = line.split()
                E_x.append(float(tmp[6]))
            except IndexError:
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
                    B_M1.append([float(tmp[1]), float(tmp[5][:-1])])

                elif (tmp[1][-1] != ")") and (tmp[3][-1] == ")"):
                    """
                    Example: 2+(10) 17.791 2+( 1) 5.172 12.619 0.006( 0.00) 0.006( 0.00)
                    """
                    B_M1.append([float(tmp[1]), float(tmp[6][:-1])])
                
                elif (tmp[1][-1] == ")") and (tmp[4][-1] != ")"):
                    """
                    Example: 3+( 8) 19.503 2+(11) 18.393 1.111 0.000( 0.00) 0.000( 0.00)
                    """
                    B_M1.append([float(tmp[2]), float(tmp[6][:-1])])

                elif (tmp[1][-1] == ")") and (tmp[4][-1] == ")"):
                    """
                    Example: 1+( 7) 19.408 2+( 9) 16.111 3.297 0.005( 0.00) 0.003( 0.00)
                    """
                    B_M1.append([float(tmp[2]), float(tmp[7][:-1])])

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
    B_M1_expected = [[5.172, 20.5], [17.791, 0.0], [19.408, 5.7], [18.393, 0.1]]
    
    for calculated, expected in zip(E_x, E_x_expected):
        msg = f"Error in E_x. Expected: {expected}, got: {calculated}."
        assert calculated == expected, msg

    for calculated, expected in zip(B_M1, B_M1_expected):
        msg = f"Error in B_M1. Expected: {expected}, got: {calculated}."
        success = (calculated[0] == expected[0]) and (calculated[1] == expected[1])
        assert success, msg



if __name__ == "__main__":
    # read_kshell_output()
    test_file_read()
    pass