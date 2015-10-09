"""Compute eigenvalues of transition matrices with varying lag times"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
import os

import tmatrix_calc as tmcalc

def calc_all(args):
    
    cwd = os.getcwd()
    
    if args.FRETdir != None:
        os.chdir(args.FRETdir)
    
    FRETfile = args.FRETfile
    FRET_trace = np.loadtxt(FRETfile)
    
    os.chdir(cwd)

    if not os.path.isdir(args.savedir):
        os.mkdir(args.savedir)

    os.chdir(args.savedir)
        
    fstep = args.fstep
    fstep = range(fstep[0], fstep[1]+fstep[2], fstep[2])
    
    if args.db:
        db = True
    else:
        db = False
    
    
    for i in range(len(fstep)):
        
        T_matrix = tmcalc.get_T_matrix(FRET_trace, spacing=args.spacing, framestep=fstep[i], flatten=False, db=db)
        np.savetxt("T_matrix_fstep_%s.dat"%fstep[i], T_matrix)
        
        lagtime = float(fstep[i])*0.5
        
        print "Calculating eigenvalues for frame step %d, lag time %1.1f ps"%(fstep[i], lagtime)
        eigs = compute_eigs(T_matrix)
        np.savetxt("Eigenvalues_fstep_%s.dat"%fstep[i], eigs, header="eigenvalues for lag time %1.1f ps"%lagtime)
        
        if args.ts:
            timescale = compute_time_scale(eigs, lagtime)
            np.savetxt("Calc_time_scales_fstep_%s.dat"%fstep[i], timescale)

    if args.tplot:
        plot_time_scale(fstep)

def compute_eigs(matrix):
    
    w,v = np.linalg.eig(matrix)
    
    #Check if there are complex eigenvalues
    num_complex = np.sum(np.iscomplex(w))
    print("Found %d complex eigenvalues."%num_complex)

    return w

def compute_time_scale(eigs, lagtime):

    if np.any(np.iscomplex(eigs)):
        raise ValueError("Not all eigenvalues are real. Cannot compute time scale.")

    #Sort eigenvalues by magnitude
    idx = np.abs(eigs).argsort()[::-1]   
    eigs = eigs[idx]
        
    #Extract useable eigenvalues
    if np.size(np.where(eigs<0)) == 0:
        eigs_use = eigs[1:]
    else:
        tr = np.min(np.where(eigs < 0))
        eigs_use = eigs[1:tr]
    
    print "The eigenvalues used for time scale calculation are: "
    print eigs_use
        
    # Compute time scale using t = -lagtime/ln(eigenvalue)
    timescale = -lagtime/np.log(eigs_use)
    
    return timescale

def plot_time_scale(fstep):

    plot_arr = np.zeros((1,1))
    
    # Construct giant matrix of values to plot
    for i in range(len(fstep)):
        tmp = np.loadtxt("Calc_time_scales_fstep_%s.dat"%fstep[i])
        
        if np.size(tmp) == 1:
            tmp = np.reshape(tmp,(np.size(tmp),))
        
        r_tmp = np.shape(tmp)[0]
        r_pta = np.shape(plot_arr)[0]
        
        # Add filler zeros if arrays have different numbers of rows
        if r_tmp < r_pta:
            tmp = np.concatenate((tmp, np.zeros(r_pta - r_tmp)))
        
        if r_tmp > r_pta:
            if np.size(np.shape(plot_arr)) == 1:
                cols = 1
            else:
                cols = np.shape(plot_arr)[1]
                
            plot_arr = np.concatenate((plot_arr, np.zeros((r_tmp - r_pta, cols))))
                
        tmp = np.reshape(tmp,(np.size(tmp),1))
        
        plot_arr = np.concatenate((plot_arr, tmp), axis=1)

    plot_arr = np.ma.masked_where(plot_arr == 0, plot_arr)

    plt.figure()
#    lagtime = 0.5*np.array(fstep, dtype=float)
    lagtime = np.array(fstep, dtype=float)
    
    for i in range(np.shape(plot_arr)[0]):
        plt.semilogy(np.transpose(lagtime), plot_arr[i,1:], 'o-', alpha=0.75, label="Eigenvalue %d"%(i+2))

    plt.axis([0, 1000, 0, 100000])
    plt.xlabel("Lag Step (frames)")
    plt.ylabel("Calculated Time Scale (frames)")
    plt.title("Time Scale vs. Lag Time")
    #plt.legend(loc=2, prop={'size':10})
    plt.savefig("time_plot.png")
    plt.close()

def get_args():

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--FRETdir", type=str, help="Location of FRET trace")
    parser.add_argument("--FRETfile", default="FRET_trace.dat", type=str, help="File containing FRET trace")
    parser.add_argument("--savedir", default="Eigenvalues", type=str, help="Save location")
    parser.add_argument("--spacing", default=0.1, type=float, help="histogram spacing")
    parser.add_argument("--fstep", nargs=3, type=int, help="Framestep range, format min, max, step. 1 frame=0.5 ps lag time")
    parser.add_argument("--db", action="store_true", help="Detailed balance matrix")
    parser.add_argument("--ts", action="store_true", help="Use if you want to compute time scale")
    parser.add_argument("--tplot", action="store_true", help="Use if you want to generate lag time vs time scale plot")
    
    args = parser.parse_args()
    
    return args

if __name__=="__main__":

    args = get_args()
    calc_all(args)