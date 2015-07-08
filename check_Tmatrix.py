import numpy as np

import argparse
import os

import tmatrix_calc as tmcalc
import plot_tmatrix as ptm
import extract_tmatrix as xtm

def plot_single_matrix(args):
    
    cwd = os.getcwd()
    subdir = args.subdir

    os.chdir(subdir)
    
    FRETfile = args.trace
    FRETr = np.loadtxt(FRETfile)
    
    # Calculate transition matrix
    tmatrix = tmcalc.get_T_matrix(FRETr, framestep=args.fstep)
    
    # Unflatten matrix if necessary
    if tmcalc.check_matrix(tmatrix):
        tmatrix = tmcalc.unflatten_matrix(tmatrix)

    # Plot matrix
    if "color" in args.plot_type:
        diff = False
        ptm.plot_matrix(tmatrix, "Transition Matrix", 0, 0.5, diff)
    if "spy" in args.plot_type:
        ptm.spy_matrix(tmatrix, "Non-Zero Transitions")

def plot_compare(args):
    
    cwd = os.getcwd()
    subdir = args.subdir
    os.chdir(subdir)
    
    FRETfile = args.trace
    FRETr = np.loadtxt(FRETfile)
    
    # Load full experimental transition matrix
    TMfile = args.TMfile
    T_matrix = np.loadtxt(TMfile)
    
    # Calculate and unflatten transition matrix
    tmatrix = tmcalc.get_T_matrix(FRETr, framestep=args.fstep)
    if tmcalc.check_matrix(tmatrix):
        tmatrix = tmcalc.unflatten_matrix(tmatrix)
    
    # Calculate experimental matrix
    bin_size = np.size(tmatrix)**0.5
    spacing=0.1
    lower_ran = spacing*(np.floor(10*np.min(FRETr)))
    upper_ran = spacing*(np.ceil(10*np.max(FRETr)))
    ran_size = (lower_ran, upper_ran)
    tmatrix_exp = xtm.tmatrix_exp_calc(T_matrix, bin_size, ran_size, spacing)
    if tmcalc.check_matrix(tmatrix_exp):
        tmatrix_exp = tmcalc.unflatten_matrix(tmatrix_exp)
    
    # Difference matrix
    diff_matrix = tmatrix - tmatrix_exp
    
    # Generate plots
    if "color" in args.plot_type:
        ptm.plot_matrix(tmatrix, "Simulated Transition Matrix", 0, 0.5)
        ptm.plot_matrix(tmatrix_exp, "Experimental Transition Matrix", 0, 0.5)
        ptm.plot_matrix(diff_matrix, "Transition Probability Differences", -0.5, 0.5, diff=True)
    if "spy" in args.plot_type:
        ptm.spy_matrix(tmatrix, "Non-zero entries, simulated")
        ptm.spy_matrix(tmatrix_exp, "Non-zero entries, experimental")
        ptm.spy_matrix(diff_matrix, "Non-zero entries, difference")

def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("--subdir", default="Tmatrix_Long", type=str, help="Subdirectory")
    parser.add_argument("--trace", default="FRET_trace.dat", type=str, help="File containing FRET trace")
    parser.add_argument("--TMfile", default="T_matrix_exp.dat", type=str, help="Experimental T matrix")
    parser.add_argument("--fstep", default=4, type=int, help="Frames between transitions, corresponds to lag time")
    parser.add_argument("--compare", action="store_true", help="Specify that you want to compare it to the experimental matrix")
    parser.add_argument("--plot_type", default="color", type=str, choices=["spy", "color"], help="Type of plot")

    args = parser.parse_args()

    return args

if __name__=="__main__":
    
    args = get_args()
    
    if args.compare:
        plot_compare(args)
    else:
        plot_single_matrix(args)