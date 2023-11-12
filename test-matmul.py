import subprocess

binary = "./matmul.exe"
matrices = [["./matrices/testMatrices/A.in", "./matrices/testMatrices/B.in", "./matrices/testMatrices/C.out"],
            ["./matrices/testMatrices/A2.in", "./matrices/testMatrices/B2.in", "./matrices/testMatrices/C2.out"],
            ["./matrices/testMatrices/A3.in", "./matrices/testMatrices/B3.in", "./matrices/testMatrices/C3.out"]]


def run_matmul(mat1, mat2, outfile):
    print(f"running: {outfile}")
    result = subprocess.run([binary, mat1, mat2, outfile, "OPT2"], capture_output=True, text=True)
    print(f"{result.stdout}")


def compare(answer, solution):
    with open(answer) as ans:
        ans_arr = ans.readlines()
    with open(solution) as sol:
        sol_arr = sol.readlines()
    assert ans_arr == sol_arr, "Not the same"


if __name__ == '__main__':
    counter = 1
    for mat in matrices:
        out = f"./out{counter}.txt"
        run_matmul(mat[0], mat[1], out)
        out = f"out{counter}.txt"
        compare(out, mat[2])
        counter += 1
