import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import kshell_utilities as ks

def trunctaion():
    data_files = []
    path = "kshell_output/run_01"
    for directory in sorted(os.listdir(path)):
        if os.path.isdir(f"{path}/{directory}"):
            data_files.append(ks.loadtxt(fname=f"{path}/{directory}"))
            
            # for summary_file in sorted(os.listdir(f"{path}/{directory}")):
            #     if summary_file.startswith("summary"):
            #         data = ks.loadtxt(fname=f"{path}/{directory}/{summary_file}")
            #         data_files.append(data)
            # for ptn_file in sorted(os.listdir(f"{path}/{directory}")):
            #     if ptn_file.endswith(".ptn"):
            #         data = ks.loadtxt(fname=f"{path}/{directory}/{ptn_file}")
            #         data_files.append(data)
            #         break

    # print(f"{data_files[1].proton_partition=}")
    # print(f"{data_files[1].neutron_partition=}")
    # print(f"{data_files[0].truncation=}")
    # path = "/Users/jon/Library/Mobile Documents/com~apple~CloudDocs/UiO/2020/FYS5960/FYS5960/exercises/truncation/src/kshell_output/run_1/trunc_01/Ni56_gxpf1a_p.ptn"
    # data = ks.loadtxt(fname=path)
    # print(f"{data.truncation=}")
    # print(f"{data.proton_partition}")
    # print(f"{data.neutron_partition}")


    line_width = 0.4
    x_scope = range(len(data_files[1].E_x))
    for i in range(len(data_files)):
        try:
            n = len(data_files[i].E_x)
        except TypeError:
            """
            E_x might be None.
            """
            continue
        
        for j in range(n):
            plt.hlines(
                y = data_files[i].E_x[j],
                xmin = x_scope[j] - line_width,
                xmax = x_scope[j] + line_width,
                color = list(mcolors.TABLEAU_COLORS.items())[i][1]
            )
        
        plt.plot([], [], label = data_files[i].truncation[0][1])    # Only for setting correct labels.
    
    plt.legend(loc="lower right")
    plt.ylabel("MeV")
    plt.show()
                    


if __name__ == "__main__":
    trunctaion()