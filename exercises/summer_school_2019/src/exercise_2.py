import numpy as np
import matplotlib.pyplot as plt
from exercise_1_levels import level_density

T = 1.238,   # [MeV].
E0 = 0.448    # [MeV].

def ct(Ex):
    return np.exp((Ex - E0)/T)/T

# if __name__ == "__main__":
#     bin_size_1 = 200    # [keV].
#     bin_size_2 = 500    # [keV].
#     counts_1, density_1, levels = level_density(bin_size_1)
#     counts_2, density_2, levels = level_density(bin_size_2)
#     Ex = np.linspace(0, 11200, 100)

#     ct_density = ct(Ex=levels)

#     # plt.plot(levels, ct_density)
#     plt.plot(levels, density_1)
#     plt.show()

