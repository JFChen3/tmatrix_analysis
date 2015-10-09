import numpy as np
import scipy.stats as stats

def find_FRET_bins(FRETr, spacing):
    # Histograms the trace data into macrostate bins, analogous to find_sim_bins from compute_Jacobian script
    
    weights = np.ones(np.shape(FRETr)[0])
    
    # Taken from find_sim_bins, included for consistency
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = 0
#    minvalue = int(np.amin(FRETr)/spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*spacing,maxvalue*spacing)

    hist, bins, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)

    return hist, slices, num_bins, ran_size, spacing

def tmatrix_exp_calc(T_matrix, bin_size, ran_size, spacing):

    #T_matrix is actually the experimental T_matrix
        
    lower_bin = int(np.around(ran_size[0]/spacing, 0))
    upper_bin = int(np.around(ran_size[1]/spacing, 0))
    
    T_matrix_small = T_matrix[lower_bin:upper_bin, lower_bin:upper_bin]
    
    print np.shape(T_matrix_small)
    np.savetxt("T_matrix_small.dat", T_matrix_small)
    
    return T_matrix_small

if __name__=="__main__":

    FRETr = np.loadtxt("FRETr.dat")
    
    hist, slices, num_bins, ran_size, spacing = find_FRET_bins(FRETr)
    
    print num_bins
    print ran_size
    print spacing
    
    TMfile = "T_matrix_exp.dat"
    T_matrix = np.loadtxt(TMfile)
    tmatrix_exp_calc(T_matrix, num_bins, ran_size, spacing)