"""Compute eigenvalues of transition matrices with varying lag times"""

import numpy as np
#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt

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
    
    if args.db:
        db = True
    else:
        db = False

    for i in range(len(fstep)):
        
        T_matrix = tmcalc.get_T_matrix(FRET_trace, framestep=fstep[i], flatten=False, db=db)
        np.savetxt("T_matrix_fstep_%s.dat"%fstep[i], T_matrix)
        
        lagtime = float(fstep[i])*0.5
        
        print "Calculating eigenvalues for frame step %d, lag time %1.1f ps"%(fstep[i], lagtime)
        eigs = compute_eigs(T_matrix)
        np.savetxt("Eigenvalues_fstep_%s.dat"%fstep[i], eigs, header="eigenvalues for lag time %1.1f ps"%lagtime)
    
def compute_eigs(matrix):
    
    w,v = np.linalg.eig(matrix)
    
    #Check if there are complex eigenvalues
    num_complex = np.sum(np.iscomplex(w))
    print("Found %d complex eigenvalues."%num_complex)

    return w

def get_args():

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--FRETdir", type=str, help="Location of FRET trace")
    parser.add_argument("--FRETfile", default="FRET_trace.dat", type=str, help="File containing FRET trace")
    parser.add_argument("--savedir", default="Eigenvalues", type=str, help="Save location")
    parser.add_argument("--fstep", nargs="+", type=int, help="Framesteps, 1 frame=0.5 ps lag time")
    parser.add_argument("--db", action="store_true", help="Detailed balance matrix")
    
    args = parser.parse_args()
    
    return args

if __name__=="__main__":

    args = get_args()
    calc_all(args)