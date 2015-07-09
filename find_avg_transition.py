"""Find average transition bin for each starting bin"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
import os

from tmatrix_calc import check_matrix, unflatten_matrix

def plot_all(args):
    
    cwd = os.getcwd()
    subdir = args.subdir
    os.chdir("%s/%s"%(cwd, subdir))
    
    if args.compare:
        T_matrix_1 = np.loadtxt(args.files[0])
        if check_matrix(T_matrix_1):
            T_matrix_1 = unflatten_matrix(T_matrix_1)
        bins_1, avg_j_1 = find_avg_transition(T_matrix_1)
        
        T_matrix_2 = np.loadtxt(args.files[1])
        if check_matrix(T_matrix_2):
            T_matrix_2 = unflatten_matrix(T_matrix_2)
        bins_2, avg_j_2 = find_avg_transition(T_matrix_2)
        labels = args.labels
        plot_compare(bins_1, avg_j_1, bins_2, avg_j_2, labels)
    else:
        T_matrix = np.loadtxt(args.files[0])
        if check_matrix(T_matrix):
            T_matrix = unflatten_matrix(T_matrix)        
        bins, avg_j = find_avg_transition(T_matrix)
        plot_single_vec(bins, avg_j)

def find_avg_transition(tmatrix):

    dim = np.shape(tmatrix)[0]

    avg_j = np.zeros(dim)
    binrange = np.reshape(range(dim),(dim,1))

    for i in range(dim):
        tmp = np.reshape(tmatrix[i,:],(1,dim))
        avg_j[i] = np.dot(tmp, binrange)

    return binrange, avg_j

def get_linfit(x, y):
    
    xmin = np.min(x)
    xmax = np.max(x)
    x = np.reshape(x,(np.size(x),))
    y = np.reshape(y,(np.size(y),))
#    x = x[10:]
#    y = y[10:]
    x = x[y != 0]
    y = y[y != 0]
    fit = np.polyfit(x,y,1)
    fit_fn = np.poly1d(fit)
    print fit_fn
    
    x_fit = np.arange(xmin, xmax, 0.1)
    y_fit = fit_fn(x_fit)
    
    return x_fit, y_fit, fit_fn

def plot_single_vec(binrange, avg_j):
    
    plt.figure()
    plt.plot(binrange, avg_j, 'bo', alpha=0.75)
    x_fit, y_fit, fit_fn = get_linfit(binrange, avg_j)
    plt.plot(x_fit, y_fit, 'b-', alpha=0.75, label="%s"%fit_fn)
    plt.xlabel("bin i")
    plt.ylabel("Average bin j")
    plt.title("Average transitions from bin i")
    plt.savefig("avg_j_bin.png")
    plt.close()

def plot_compare(bins_1, avg_j_1, bins_2, avg_j_2, labels):
    
    plt.figure()
    x_fit_1, y_fit_1, fit_fn_1 = get_linfit(bins_1, avg_j_1)
    plt.plot(bins_1, avg_j_1, 'bo', alpha=0.75, label="%s"%labels[0])
    plt.plot(x_fit_1, y_fit_1, 'b-', alpha=0.75, label="%s"%fit_fn_1)
    
    plt.plot(bins_2, avg_j_2, 'ro', alpha=0.75, label="%s"%labels[1])
    x_fit_2, y_fit_2, fit_fn_2 = get_linfit(bins_2, avg_j_2)
    plt.plot(x_fit_2, y_fit_2, 'r-', alpha=0.75, label="%s"%fit_fn_2)
    
    plt.xlabel("bin i")
    plt.ylabel("Average bin j")
    plt.title("Average transitions from bin i")
    plt.legend(loc=2)
    plt.savefig("avg_j_bin_compare.png")
    plt.close()

def get_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--subdir", type=str, help="Subdirectory containing files")
    parser.add_argument("--files", nargs="+", type=str, help="Files containing matrices")
    parser.add_argument("--labels", nargs="+", type=str, help="Legend labels")
    parser.add_argument("--compare", action="store_true", help="Specify that you want to compare")
    
    args = parser.parse_args()
    
    return args

if __name__=="__main__":
    
    args = get_args()
    plot_all(args)        