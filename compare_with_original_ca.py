import sys
import subprocess as sp
import re

import numpy as np

import coupled_alpha


def run_orig_ca(xs, ys):
    ca = coupled_alpha.CoupledAlpha(xs, ys)
    
    values = {}
    for k, k_lists in enumerate(ca.ST.lists):
        for key in k_lists:
            for node in k_lists[key]:
                values[tuple(ca.ST.node2simplex(node))] = np.sqrt(float(node.filtration))
    
    return values


def run_cpp_ca(xs, ys):
    dim = xs.shape[1]

    with open("labels_pc.txt", "w") as f:
        for k in range(xs.shape[0]):
            print("0 {}".format(" ".join(map(str, xs[k, :]))), file=f)

        for k in range(ys.shape[0]):
            print("1 {}".format(" ".join(map(str, ys[k, :]))), file=f)

    result = sp.run(["../coupled_alpha", str(dim), "labels_pc.txt"], capture_output=True, text=True, check=True)

    values = {}
    for line in re.split(r"\n", result.stdout):
        if line == "":
            continue
        l = re.split(r"\s+", line)
        value = l.pop()
        values[tuple(map(int, l))] = float(value)

    return values
    

def main(dim, m, n):
    xs = np.random.uniform(0, 1, (m, dim))
    ys = np.random.uniform(0, 1, (n, dim))
    values_1 = run_orig_ca(xs, ys)
    values_2 = run_cpp_ca(xs, ys)
    print(values_1.keys() == values_2.keys())
    for k in values_1:
        if abs(values_1[k] - values_2[k])  > 1e-8:
            print(k, values_1[k], values_2[k])
    
if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
