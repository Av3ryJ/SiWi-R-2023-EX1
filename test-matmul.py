import subprocess

binary = "matmul.exe"
matrices = [["./matrices/testMatrices/A.in", "./matrices/testMatrices/B.in", "./matrices/testMatrices/C.out"],
            ["./matrices/testMatrices/A1.in", "./matrices/testMatrices/B1.in", "./matrices/testMatrices/C1.out"],
            ["./matrices/testMatrices/A2.in", "./matrices/testMatrices/B2.in", "./matrices/testMatrices/C2.out"]]


def run_matmul(mat1, mat2, out):
    subprocess.run([binary, mat1, mat2, f">{out}"])


def compare(answer, solution):
    with open(answer) as ans:
        ans_arr = ans.readlines()
    with open(solution) as sol:
        sol_arr = sol.readlines()
    assert ans_arr == sol_arr, "Not the same"


if __name__ == '__main__':
    for mat in matrices:
        out = "out"
        run_matmul(mat[0], mat[1], out)
        compare(out, mat[2])
