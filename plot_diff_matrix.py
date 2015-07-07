"""Plot difference in transition matrices between two iterations"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
import os

from tmatrix_calc import check_matrix, unflatten_matrix
import plot_tmatrix as ptm

def plot_difference_matrix(args):

    subdir = args.subdir

    cwd = os.getcwd()

    os.chdir("%s/%s"%(cwd,subdir[0]))
    T_matrix_1 = np.loadtxt(args.files[0])
    if check_matrix(T_matrix_1):
        T_matrix_1 = unflatten_matrix(T_matrix_1)

    os.chdir(cwd)

    os.chdir("%s/%s"%(cwd,subdir[1]))
    T_matrix_2 = np.loadtxt(args.files[1])
    if check_matrix(T_matrix_2):
        T_matrix_2 = unflatten_matrix(T_matrix_2)

    os.chdir(cwd)

    c_direc = "%s/tmatrix_compare"%(cwd)
    
    if not os.path.isdir(c_direc):
        os.mkdir(c_direc)
    
    os.chdir(c_direc)
    
    diff_matrix = T_matrix_1 - T_matrix_2
    
    if args.iter != None:
        iter = args.iter
        title = "Transition Differences Iter %d and %d"%(iter[0], iter[1])
    else:
        title = "Transition Differences"
    
    ptm.plot_matrix(diff_matrix, title, -0.5, 0.5, diff=True)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--subdir", nargs=2, type=str, help="Directory containing files")
    parser.add_argument("--iter", nargs=2, type=int, help="Iterations to compare")
    parser.add_argument("--files", default=("T_matrix_sim.dat", "T_matrix_sim.dat"), nargs=2, type=str, help="T_matrix file names")

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = get_args()
    plot_difference_matrix(args)