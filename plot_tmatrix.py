"""Plot transition matrices"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse

from tmatrix_calc import check_matrix, unflatten_matrix

def plot_single(args):
    # Plot a single matrix

    matrix_file = args.tmatrix

    T_matrix = np.loadtxt(matrix_file)

    if check_matrix(T_matrix):
        T_matrix = unflatten_matrix(T_matrix)

    if "color" in args.plot_type:
        plot_matrix(T_matrix, "Transition Matrix", 0, 0.5)
    if "spy" in args.plot_type:
        spy_matrix(T_matrix, "Non-Zero Transitions")

    return

def plot_double(args):
    # Plot two matrices along with a difference matrix
    
    if args.iter != None:
        import os
        subdir = args.subdir
        iter = args.iter
        file_location = "%s/iteration_%d"%(subdir, iter)
        os.chdir(file_location)
        print os.getcwd()

    ##ASSUMES SIMULATED MATRIX FED FIRST
    (sim_file, exp_file) = args.tmatrix

    sim_matrix = np.loadtxt(sim_file)
    exp_matrix = np.loadtxt(exp_file)
    
    if check_matrix(sim_matrix):
        sim_matrix = unflatten_matrix(sim_matrix)
    if check_matrix(exp_matrix):
        exp_matrix = unflatten_matrix(exp_matrix)
    
    # Subtract matrices to determine difference matrix
    diff_matrix = sim_matrix - exp_matrix
    
    ##PRINT AVERAGE DIFFERENCE OF NON-ZERO ENTRIES
    print np.mean(diff_matrix)
    
    if "color" in args.plot_type:
        plot_matrix(sim_matrix, "Simulated Transition Matrix", 0, 0.5)
        plot_matrix(exp_matrix, "Experimental Transition Matrix", 0, 0.5)
        plot_matrix(diff_matrix, "Transition Probability Differences", -0.5, 0.5, diff=True)
    if "spy" in args.plot_type:
        spy_matrix(sim_matrix, "Non-zero entries, simulated")
        spy_matrix(exp_matrix, "Non-zero entries, experimental")
        spy_matrix(diff_matrix, "Non-zero entries, difference")
    
    return

def spy_matrix(matrix, title):
    # Plot non-zero matrix entries

    zmin = np.min(matrix)
    zmax = np.max(matrix)

    plt.figure()

    plt.spy(matrix, precision=0.01)

    plt.xlabel("j")
    plt.ylabel("i")
    plt.title("%s"%title)

    plt.savefig("%s"%title)
    plt.close()

def plot_matrix(matrix, title, zmin, zmax, diff=False):
    # Plot color map

    plt.figure()
    if diff == True:
        plt.pcolormesh(matrix, cmap="RdBu", vmin=zmin, vmax=zmax)
    else:
        plt.pcolormesh(matrix, vmin=zmin, vmax=zmax)
    plt.colorbar()

    plt.axis([0, np.shape(matrix)[1], 0, np.shape(matrix)[0]])
    plt.xlabel("j")
    plt.ylabel("i")
    plt.title("%s"%title)

    plt.savefig("%s"%title)
    plt.close()

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--tmatrix", default=("sim_feature.dat", "target_feature.dat"), nargs='+', type=str, help="File containing T-matrix")
    parser.add_argument("--plot_type", default="color", nargs='+', type=str, choices=["spy", "color"], help="Type of plot")
    parser.add_argument("--subdir", type=str, required=False, help="Subdirectory")
    parser.add_argument("--iter", type=int, required=False, help="Iteration to plot")

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = get_args()

    matrix_file = args.tmatrix

    if len(matrix_file) == 1:
        plot_single(args)
    elif len(matrix_file) == 2:
        plot_double(args)
    else:
        print "ERROR: Too many input files!"