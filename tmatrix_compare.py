import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tmatrix_analysis.tmatrix_calc as tmcalc

import os

def track_over_iterations(iter_min, iter_max):
    
    iter_range = range(iter_min, iter_max+1)
    
    cwd = os.getcwd()
    
    num_iters = iter_max - iter_min + 1
    abs_diff = np.zeros(num_iters)
    sq_diff = np.zeros(num_iters)
    dot_product = np.zeros(num_iters)
    
    for i in range(num_iters):
        os.chdir("%s/iteration_%d"%(cwd, iter_range[i]))
        matrix_1 = np.loadtxt("sim_feature.dat")
        matrix_2 = np.loadtxt("target_feature.dat")
        diff_matrix = calc_diff_matrix(matrix_1, matrix_2, nonzero=True)
        abs_diff[i] = avg_abs_diff(diff_matrix)
        sq_diff[i] = avg_sq_diff(diff_matrix)
        dot_product[i] = vector_dot_product(matrix_1, matrix_2)
    
    print "\nSim Feature vs. Target Feature Comparisons, Iterations %d to %d"%(iter_min, iter_max)
    print "\nAverage absolute difference: "
    print abs_diff
    print "\nAverage squared difference: "
    print sq_diff
    print "\nDot product: "
    print dot_product
    
    os.chdir(cwd)
    
    if not os.path.isdir("%s/sim_exp_compare"%cwd):
        os.mkdir("%s/sim_exp_compare"%cwd)

    os.chdir("%s/sim_exp_compare"%cwd)
    
    # Save data to a file
    file = "sim_exp_compare_iter_%d_%d.txt"%(iter_min, iter_max)
    f = open(file,"w")
    f.write("Sim Feature vs. Target Feature Comparisons, Iterations %d to %d"%(iter_min, iter_max))
    f.write("\n")
    f.write("\nAverage absolute difference: \n")
    f.write(np.array_str(abs_diff))
    f.write("\nAverage squared difference: \n")
    f.write(np.array_str(sq_diff))
    f.write("\nDot product: \n")
    f.write(np.array_str(dot_product))
    f.close()
    
    plt.figure()
    plt.plot(range(iter_min, iter_max+1), abs_diff)
    
def calc_diff_matrix(matrix_1, matrix_2, nonzero=True):
    
    diff_matrix = matrix_1 - matrix_2
    if nonzero:
        diff_masked = np.ma.masked_where(matrix_1+matrix_2==0, diff_matrix)
        diff_matrix = diff_masked
    
    return diff_matrix

def avg_abs_diff(diff_matrix):
    # Find average square of a difference matrix
    
    abs_diff_matrix = np.abs(diff_matrix)
    abs_diff = np.mean(abs_diff_matrix)
    
    return abs_diff

def avg_sq_diff(diff_matrix):

    sq_matrix = np.square(diff_matrix)
    sq_diff = np.mean(sq_matrix)
    
    return sq_diff

def vector_dot_product(matrix_1, matrix_2):
    
    if not tmcalc.check_matrix(matrix_1):
        matrix_1 = np.ndarray.flatten(matrix_1)
    matrix_1 /= np.linalg.norm(matrix_1)
    if not tmcalc.check_matrix(matrix_2):
        matrix_2 = np.ndarray.flatten(matrix_2)
    matrix_2 /= np.linalg.norm(matrix_2)
    
    dot_product = np.dot(matrix_1, matrix_2)
    
    return dot_product

if __name__=="__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--subdir", type=str, help="Directory containing iteration_N folders")
    parser.add_argument("--iter", nargs=2, type=int, help="Iterations to track, in format [iter_min iter_max]")
    
    args = parser.parse_args()
    
    cwd = os.getcwd()
    
    os.chdir("%s"%args.subdir)
    
    iter_min = args.iter[0]
    iter_max = args.iter[1]
    
    track_over_iterations(iter_min, iter_max)