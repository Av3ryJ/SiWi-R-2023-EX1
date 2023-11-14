import subprocess
import os.path as p
import json
import matplotlib.pyplot as plt
import numpy as np

binary = "./matmul.exe"
sizes = [32, 64, 128, 256, 512, 1024, 2048]
options = ["STD", "BLAS", "OPT1", "OPT2"]
matrix_folder = "./matrices/perfMatrices"
json_name = "./times.json"
loaded_json = {option: {size: 0. for size in sizes} for option in options}


def run_matmul(mat1, mat2, outfile, option):
    result = subprocess.run([binary, mat1, mat2, outfile, option], capture_output=True, text=True)
    return result.stdout


def time_all():
    for size in sizes:
        for option in options:
            mat = matrix_folder + f"/{size}x{size}-"
            returned = run_matmul(mat+"1", mat+"2", "out.txt", option)
            time = returned.split()
            loaded_json[option][size] = float(time[0])
    with open(json_name, 'w') as f:
        json.dump(loaded_json, f)


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

    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(np.log2(sizes), get_array_for_option("STD"))
    axs[0, 0].set_ylim(bottom=0, top=max_y)
    axs[0, 0].set_title("Option: STD")

    axs[0, 1].plot(np.log2(sizes), get_array_for_option("BLAS"))
    axs[0, 1].set_ylim(bottom=0, top=max_y)
    axs[0, 1].set_title("Option: Blas")

    axs[1, 0].plot(np.log2(sizes), get_array_for_option("OPT1"))
    axs[1, 0].set_ylim(bottom=0, top=max_y)
    axs[1, 0].set_title("Option: OPT1")

    axs[1, 1].plot(np.log2(sizes), get_array_for_option("OPT2"))
    axs[1, 1].set_ylim(bottom=0, top=max_y)
    axs[1, 1].set_title("Option: OPT2")

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    if p.exists(json_name):
        with open(json_name, 'r') as file:
            loaded_json = json.load(file)
    else:
        time_all()

    plot_by_option()


