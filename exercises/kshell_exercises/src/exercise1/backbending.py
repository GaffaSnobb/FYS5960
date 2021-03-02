import numpy as np
import matplotlib.pyplot as plt

def read_file():
    fname = "summary_Cr48_gxpf1a.txt"
    J = []
    E = []

    with open(fname, "r") as infile:
        for _ in range(5): infile.readline()
        for line in infile:
            try:
                tmp = line.split()
                J.append(float(tmp[1]))
                E.append(float(tmp[6]))

            except IndexError:
                continue


    J = np.array(J)
    E = np.array(E)
    print(E)

    plt.plot(E[1:] - E[:-1], J[1:])
    plt.show()


read_file()