import subprocess
import os.path as p
import json
import matplotlib.pyplot as plt
import numpy as np

# stuff for plots
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
arrowprops=dict(arrowstyle="->", connectionstyle="arc")
kw = dict(xycoords='data', textcoords="axes fraction",
          arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")


blocksizes = [16, 32, 64, 128, 256, 512, 1024]
json_for_strassen = "./strassen.json"
loaded_strassen = {"OPT2": [0. for j in blocksizes], "OPT3": [0. for i in blocksizes]}

binary = "./matmul.exe"
sizes = [32, 64, 128, 256, 512, 1024, 2048]
options = ["STD", "BLAS", "OPT1", "OPT2", "OPT3"]
matrix_folder = "./matrices/perfMatrices"
json_name = "./times.json"
loaded_json = {option: {size: 0. for size in sizes} for option in options}


def run_matmul(bin, mat1, mat2, outfile, option):
    result = subprocess.run([bin, mat1, mat2, outfile, option], capture_output=True, text=True)
    return result.stdout


def time_all():
    for size in sizes:
        for option in options:
            mat = matrix_folder + f"/{size}x{size}-"
            returned = run_matmul(binary, mat+"1", mat+"2", "out1.txt", option)
            time = returned.split()
            loaded_json[option][size] = float(time[0])
    with open(json_name, 'w') as f:
        json.dump(loaded_json, f)


def time_strassen(size: int):
    for option in ["OPT2", "OPT3"]:
        counter = 0
        for blocksize in blocksizes:
            binn = f"./matmul{blocksize}.exe"
            mat = matrix_folder + f"/{size}x{size}-"
            returned = run_matmul(binn, mat+"1", mat+"2", "out1.txt", option)
            time = returned.split()
            loaded_strassen[option][counter] = float(time[0])
            counter += 1
    with open(json_for_strassen, 'w') as f:
        json.dump(loaded_strassen, f)


def get_array_for_option(option):
    out = []
    for key, value in loaded_json[option].items():
        out.append(value)
    return out


def plot_by_option():
    max_y = 0
    for option in options:
        max_y = max(max_y, max(get_array_for_option(option)))
    max_y *= 1.1

    fig, axs = plt.subplots(2, 1)

    array = get_array_for_option("STD")
    axs[0].plot(np.log2(sizes), array, label="STD")
    axs[0].set_ylim(bottom=0, top=max_y)
    axs[0].set_title("Runtime for different Implementations and matrix sizes")

    array = get_array_for_option("BLAS")
    axs[0].plot(np.log2(sizes), array, label="BLAS")

    array = get_array_for_option("OPT1")
    axs[0].plot(np.log2(sizes), array, label="Transpose")

    array = get_array_for_option("OPT2")
    axs[0].plot(np.log2(sizes), array, label="Strassen")

    array = get_array_for_option("OPT3")
    axs[0].plot(np.log2(sizes), array, label="Strassen + Transposed")

    array = loaded_strassen["OPT2"]
    axs[1].plot(np.log2(blocksizes), array, label="Strassen + STD")
    array = loaded_strassen["OPT3"]
    axs[1].plot(np.log2(blocksizes), array, label="Strassen + Transposed")
    axs[1].set_ylim(bottom=12, top=30)
    axs[1].set_title("Speed of Strassen with different minimal sizes")
    axs[1].annotate(f"min at 5^2 matrices", xy=(5, array[1]), xytext=(0.3, 0.5), **kw)

    for i in range(2):
        axs[i].set_xlabel("Matrix size in log2")
        axs[i].legend(loc="upper left")

    axs[0].set_ylabel("Runtime in seconds")


    #axs[2, 1].set_xlabel("Smallest size for division in log2")

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    if p.exists(json_name):
        with open(json_name, 'r') as file:
            loaded_json = json.load(file)
    else:
        time_all()
    if p.exists(json_for_strassen):
        with open(json_for_strassen, 'r') as file:
            loaded_strassen = json.load(file)
    else:
        time_strassen(2048)

    plot_by_option()


