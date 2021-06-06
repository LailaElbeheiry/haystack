import subprocess
import os
import re
# import matplotlib
# import matplotlib.pyplot as plt
# import numpy as np
# import json

polybench_dir = "/home/laila/Desktop/projects/thesis/FemtoClouds/evaluation/polybench-c-4.2.1-beta"
utilities_dir = "{}/utilities".format(polybench_dir)
polybench_spec = "{}/polybench.spec".format(utilities_dir)

dataset_size = 3
cache_size = 1024
line_size = 8
assoc = 8

poly_bins = {
    'correlation'    : '{}/datamining/correlation'.format(polybench_dir),
    'covariance'     : '{}/datamining/covariance'.format(polybench_dir),
    '3mm'            : '{}/linear-algebra/kernels/3mm'.format(polybench_dir),
    '2mm'            : '{}/linear-algebra/kernels/2mm'.format(polybench_dir),
    'atax'           : '{}/linear-algebra/kernels/atax'.format(polybench_dir),
    'bicg'           : '{}/linear-algebra/kernels/bicg'.format(polybench_dir),
    'doitgen'        : '{}/linear-algebra/kernels/doitgen'.format(polybench_dir),
    'mvt'            : '{}/linear-algebra/kernels/mvt'.format(polybench_dir),
    'gemm'           : '{}/linear-algebra/blas/gemm'.format(polybench_dir),
    'gemver'         : '{}/linear-algebra/blas/gemver'.format(polybench_dir),
    'gesummv'        : '{}/linear-algebra/blas/gesummv'.format(polybench_dir),
    'symm'           : '{}/linear-algebra/blas/symm'.format(polybench_dir),
    'syr2k'          : '{}/linear-algebra/blas/syr2k'.format(polybench_dir),
    'syrk'           : '{}/linear-algebra/blas/syrk'.format(polybench_dir),
    'trmm'           : '{}/linear-algebra/blas/trmm'.format(polybench_dir),
    'cholesky'       : '{}/linear-algebra/solvers/cholesky'.format(polybench_dir),
    'durbin'         : '{}/linear-algebra/solvers/durbin'.format(polybench_dir),
    'gramschmidt'    : '{}/linear-algebra/solvers/gramschmidt'.format(polybench_dir),
    'lu'             : '{}/linear-algebra/solvers/lu'.format(polybench_dir),
    'ludcmp'         : '{}/linear-algebra/solvers/ludcmp'.format(polybench_dir),
    'trisolv'        : '{}/linear-algebra/solvers/trisolv'.format(polybench_dir),
    'deriche'        : '{}/medley/deriche'.format(polybench_dir),
    # 'floyd-warshall' : '{}/medley/floyd-warshall'.format(polybench_dir),
    # 'nussinov'       : '{}/medley/nussinov'.format(polybench_dir),
    'adi'            : '{}/stencils/adi'.format(polybench_dir),
    'fdtd-2d'         : '{}/stencils/fdtd-2d'.format(polybench_dir),
    'heat-3d'        : '{}/stencils/heat-3d'.format(polybench_dir),
    'jacobi-1d'      : '{}/stencils/jacobi-1d'.format(polybench_dir),
    'jacobi-2d'      : '{}/stencils/jacobi-2d'.format(polybench_dir),
    'seidel-2d'      : '{}/stencils/seidel-2d'.format(polybench_dir)
}

results={}
sizes={}

def parse_output(output):
    output = output.replace(",", "")
    pat = r"\d+\.?\d*"
    res = []
    for line in output.strip().split("\n"):
        res.append(re.findall(pat, line)[0])
    return [float(i) for i in res]

def params(app):
    output = ""
    s = sizes[app]
    for param in s:
        output += " {}={}".format(param, s[param])
    return output

with open(polybench_spec) as f:
    lines = f.readlines()
    for line in lines[1:]:
        res = {}
        line = line.strip().split()
        num_vars = (len(line) - 3) // 6
        for i in range(num_vars):
            var_name = line[3 + i].lower()
            var_val = int(line[3 + num_vars + dataset_size * num_vars + i])
            res[var_name] = var_val
        sizes[line[0]] = res
    f.close()

for b in poly_bins:
    cmd = "./haystack -f {}/{}.c -c {} -l {} -k {} -I {} -d {} | grep \"TIME_TAKEN\\|COMPULSORY\\|CAPACITY\" ".format(poly_bins[b], b, cache_size, line_size, assoc, utilities_dir, params(b))
    res = subprocess.check_output(cmd, shell=True)
    res = parse_output(res.decode("utf-8"))
    res[1] = res[1] + res.pop()
    print(b, res)

    results[b] = res

with open('results.txt', 'w') as f:
    for b in poly_bins:
        f.write(b)
        f.write(',')
    f.write('\n')
    for b in poly_bins:
        f.write(str(results[b][0]))
        f.write(',')
    f.write('\n')
    for b in poly_bins:
        f.write(str(results[b][1]))
        f.write(',')
    f.write('\n')
    f.close()



