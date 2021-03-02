import numpy as np
import matplotlib.pyplot as plt

def level_density(bin_size):

    levels = np.loadtxt("levels_56Fe.txt")
    n_levels = len(levels)
    counts = np.zeros(n_levels)

    for i in range(n_levels):
        for j in range(n_levels):
            if (levels[j] >= levels[i] - bin_size) and (levels[j] <= levels[i] + bin_size):
                counts[i] += 1

    return counts, counts/bin_size, levels

# if __name__ == "__main__":
#     bin_size_1 = 200    # [KeV].
#     bin_size_2 = 500    # [KeV].
#     counts_1, density_1, levels = level_density(bin_size_1)
#     counts_2, density_2, levels = level_density(bin_size_2)

#     # plt.plot(levels, density_1)
#     # plt.plot(levels, density_2)
#     # plt.bar(levels, counts_1, width=50)
#     plt.step(levels, (density_1))
#     plt.step(levels, (density_2))
#     plt.show()