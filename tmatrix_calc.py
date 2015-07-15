"""
Calculate a transition matrix

"""

import numpy as np
import scipy.stats as stats

def find_FRET_bins(FRETr, fstep):
    # Histograms the trace data into macrostate bins
    
    weights = np.ones(np.shape(FRETr)[0])
    spacing = 0.1
    
    # Taken from find_sim_bins, included for consistency
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = 0
#    minvalue = int(np.amin(FRETr)/spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*spacing,maxvalue*spacing)
    
    hist, bins, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
    
    # Call get_transition_bins to make transition bin indices
    F_indices, t_indices = get_transition_bins(slices, num_bins, framestep=fstep)
    
    return hist, F_indices, t_indices, num_bins
    
def get_transition_bins(slices, num_bins, framestep=4):
    # Returns linear indices of "transition bins," should match with numpy.ndarray.flatten() results
    # 'slices' has dimension Nx1, where N=number of frames
    
    # Get i and j slices, i=macrostate before transition, j=macrostate after transition
    # Subtract one to convert to zero-based indices
    F_indices = slices - 1
    state_i = F_indices[:-framestep]
    state_j = F_indices[framestep:]
    
    # Compute bin indices for transitions
    t_indices = state_i*num_bins + state_j
    
    print "Bin indices are: "
    print F_indices
    print t_indices
    
    return F_indices, t_indices
    
def get_T_matrix(FRET_trace, framestep=4, flatten=False, db=False, sliding=True):
    # Calculate flattened transition matrix for a given FRET trace
    # Based on fret_analysis/compute_transitions.py
    
    # Get FRET bins
    hist, F_indices, t_indices, num_bins = find_FRET_bins(FRET_trace, framestep)
    
    T_matrix = np.zeros((num_bins, num_bins))
    
    if sliding:
        ran = range(np.shape(F_indices)[0] - framestep)
    else:
        ran = range(0, np.shape(F_indices)[0] - framestep, framestep)
    
    # Add ones to transition bins in square transition matrix
    for i in ran:
        T_matrix[F_indices[i], F_indices[i+framestep]] += 1
    
    np.savetxt("count_matrix_fstep_%d.dat"%framestep, T_matrix)
    if db:
        # Calculate detailed balance matrix
        import pyemma
        T_matrix = trim_matrix(T_matrix)
        T_matrix = pyemma.msm.estimation.transition_matrix(T_matrix, reversible=True)
    else:
        # Mask zeros to avoid divide-by-zero in normalization
        T_masked = np.ma.masked_where(T_matrix == 0, T_matrix)

        # Normalize each row
        for i in range(np.shape(T_matrix)[0]):
            T_masked[i,:] /= np.sum(T_masked[i,:])
            T_matrix = T_masked.filled(0)
    
    np.savetxt("T_matrix_sim.dat", T_matrix)
    
    # Reshape to column vector
    if flatten:
        T_matrix = np.transpose(np.ndarray.flatten(T_matrix))

    return T_matrix

def trim_matrix(matrix):
    # Remove states to avoid connectivity issues when calculating detailed balance matrix
    
    empty_bins = np.where(~matrix.any(axis=1))[0]
    print "Found %d empty bins"%np.size(empty_bins)
    
    num_states = np.shape(matrix)[0]
    
    s_lower = empty_bins[empty_bins < int(num_states/2)]
    s_upper = empty_bins[empty_bins > int(num_states/2)]
    
    if np.size(s_upper) != 0:
        tr = np.min(s_upper)
        matrix = matrix[:tr,:tr]
        print "Removed %d upper bins"%(num_states - tr)
    
    if np.size(s_lower) != 0:
        tr = np.max(s_lower) + 1
        matrix = matrix[tr:,tr:]
        print "Removed %d lower bins"%(tr)
    
    return matrix

def check_matrix(matrix):
    # Check if matrix needs to be unflattened
    
    unflatten = np.shape(matrix)[0] == np.size(matrix)
    
    return unflatten

def unflatten_matrix(matrix):
    # Unflatten matrix

    dim = np.shape(matrix)[0]**0.5

    T_matrix = np.zeros((dim, dim))

    matrix = np.transpose(matrix)
    for i in range(int(dim)):
        T_matrix[i,:] = matrix[dim*i:dim*(i+1)]

    return T_matrix

if __name__=="__main__":

    trace = np.loadtxt("FRETr.dat")
    print np.min(trace)
    print np.max(trace)
    
    T_vec = get_T_matrix(trace)
    print np.shape(T_vec)
    np.savetxt("T_vector.dat", T_vec)
            
    hist, F_indices, t_indices, num_bins = find_FRET_bins(trace)