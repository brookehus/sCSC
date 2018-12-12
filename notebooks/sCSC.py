# Python 3 file to generate sCSC dendrogram
# B.E. Husic, K.L. Schlueter-Kuck, and J.O. Dabiri
# Stanford University
# 2017-2018
# Instructions: Execute this code in the terminal, e.g.
# >>> python3 sCSC.py protein_g_adjmat.npy 7 my_directory
# Adjacency matrix must be in numpy format

import numpy as np
import matplotlib.pyplot as plt
import os, sys, pickle

# The following was written using SciPy Version 0.19.1
import scipy.linalg
from scipy.cluster.hierarchy import average, linkage, fcluster
from scipy.spatial.distance import pdist

# Compute the color field and obtains the
# user-specified number of eigenvectors
def get_eigenvectors_for_dendrogram(adjacency_matrix,
                                    numeigvecs):
    A = adjacency_matrix

    # Create column vector of row-sums of adjacency matrix
    Adegree = np.sum(A, axis=1)

    # Create diagonal degree matrix with row-sums of adjacency matrix
    DD = np.diag(Adegree)

    # Create graph Laplacian
    L = DD - A 

    # Compute eigenvectors and eigenvalues of generalized eigenvalue problem. 
    eigval, eigvec = scipy.linalg.eig(L, DD)

    # Sort eigenvalues in decending order
    lambdaindex = np.argsort(eigval)[::-1]
    sortedeigs = eigval[lambdaindex]

    # Sort eigenvectors in order corresponding to sorted eigenvalues
    CSC_fullset = eigvec[:, lambdaindex]
    
    # Select number of eigenvectors to include in dendrogram analysis
    evec = CSC_fullset[:,:numeigvecs]
    
    # This normalization makes results consistent with matlab
    for ind, ev in enumerate(evec.T):
        constant = np.matmul(np.matmul(ev.T, DD), ev)
        evec[:,ind] = ev/np.sqrt(constant)
    
    return evec, L

# Obtain binary codes for each orthogonal process
# using agglomerative hierarchical clustering
def get_binary_codes(evec, numeigvecs):
    binary = np.empty(evec.shape)

    #Loop through eigenvectors to create binary code for each state in dataset
    for col in np.arange(numeigvecs):

        # Create agglomerative hierarchical tree using average linkage
        p = pdist([[x] for x in evec[:,col]])
        lnk = linkage(p, 'average')

        # Assign value of 1 or 2 to each state state depending on membership
        # in bifurcated tree
        T = fcluster(lnk, criterion='maxclust', t=2) - 1

        binary[:, col] = T
        
    return binary

# Initialize the matrices that will be used
# to plot the dendrogram and a dictionary that will
# store cluster assignments
def initialize_dendrogram_matrices(numeigvecs):
    # Initialize matrix containing length of each dendrogram branch
    groupdist = np.zeros([2 ** (numeigvecs),numeigvecs])

    # Initialize matrix containing x coordinates of each dendrogram node
    plotcoordsx = np.zeros([2 ** (numeigvecs),numeigvecs])

    # Initialize matrix containing y coordinates of each dendrogram node
    plotcoordsy = np.zeros([2 ** (numeigvecs),numeigvecs])
    
    return groupdist, plotcoordsx, plotcoordsy, {}

# Get subgroup information
def get_subgroup_similarity(evec, L, binary, combos, split, index_dict):
    n = evec.shape[0]
    
    # Create binary code for group 1 of a pair of branches
    group1 = np.array([int(x) for x in np.binary_repr(combos, split+1)])

    # Create binary code for group 2 of a pair of branches
    group2 = np.array([int(x) for x in np.binary_repr(combos+1, split+1)])

    # Create string-compatible version of binary code for group 1 for plotting
    g1label = np.binary_repr(combos, split+1)

    # Create strong-compatible version of binary code for group 2 for plotting
    g2label = np.binary_repr(combos+1, split+1)

    # Initialize eigenvector for current level of dendrogram
    Xgroup = evec[:,split]

    # Remove eigenvector elements not in in group 1 or group 2
    not_in_g1 = np.delete(np.arange(n), np.where((binary[:,:split+1][:]==group1).all(1)))
    not_in_g2 = np.delete(np.arange(n), np.where((binary[:,:split+1][:]==group2).all(1)))
    Xgroup = np.delete(Xgroup, np.intersect1d(not_in_g1, not_in_g2))

    # Initialize graph Laplacian 
    Lgroup = L

    # Remove columns not in group 1 or group 2
    Lgroup = np.delete(Lgroup, np.intersect1d(not_in_g1, not_in_g2), axis=0)

    # Remove rows not in group 1 or group 2
    Lgroup = np.delete(Lgroup, np.intersect1d(not_in_g1, not_in_g2), axis=1)

    # Calculate number of states in group 1
    Z1count = sum((binary[:,:split+1][:]==group1).all(1))

    # If group 1 is occupied by any states
    if Z1count > 0:
        # Calculate dissimilarity metric for states in groups 1 and 2
        Z1 = np.matmul(np.matmul(Xgroup.transpose(), Lgroup), Xgroup)
        
        index_dict[g1label] = [i for i, x in enumerate(
                               (binary[:,:split+1][:]==group1).all(1))
                               if x]

    else:
        # dissimilarity metric is undefined if group 1 is unoccupied
        Z1 = np.nan

    # Calculate number of states in group 2
    Z2count = sum((binary[:,:split+1][:]==group2).all(1))

    # if group 2 is occupied by any states
    if Z2count > 0:
        # Calculate dissimilarity metric for states in groups 1 and 2
        Z2 = np.matmul(np.matmul(Xgroup.transpose(), Lgroup), Xgroup)
        
        index_dict[g2label] = [i for i, x in enumerate(
                               (binary[:,:split+1][:]==group2).all(1))
                               if x]

    else:
        # Dissimilarity metric is undefined if group 2 is unoccupied
        Z2 = np.nan
    
    return g1label, g2label, Z1count, Z2count, Z1, Z2, index_dict

# Helper function for plotting lines
def plot_line(x1, x2, y1, y2, count, adjacency_matrix, maxlinewidth,
              style='-o', color='black', markerfacecolor='black',
              markeredgecolor='black'):
    plt.plot([x1, x2], [y1, y2], style, color=color,
             markerfacecolor=markerfacecolor,
             markeredgecolor=markeredgecolor,
             linewidth=maxlinewidth*(count/len(adjacency_matrix[:,1])))

# Plot the next level of the dendrogram
def plot_next_dendrogram_split(plotcoordsx, plotcoordsy, row, split, adjacency_matrix,
                               g1label, g2label, Z1count, Z2count, Z1, Z2,
                               figsize=(14,8), lineangle=np.pi/4, maxlinewidth=15):
    
    A = adjacency_matrix
    
    xrel1 = -Z1*np.sin(lineangle)
    xrel2 = Z2*np.sin(lineangle)

    yrel1 = -Z1*np.cos(lineangle)
    yrel2 = -Z2*np.cos(lineangle)
    
    if split == 0:
        plotcoordsx[row, split] = xrel1
        plotcoordsx[row+1, split] = xrel2

        plotcoordsy[row, split] = yrel1
        plotcoordsy[row+1, split] = yrel2

        fig, ax = plt.subplots(figsize=figsize)

        plot_line(0, plotcoordsx[row, split], 0, plotcoordsy[row, split],
                  Z1count, adjacency_matrix, maxlinewidth)
        
        plot_line(0, plotcoordsx[row+1, split], 0, plotcoordsy[row+1, split],
                  Z2count, adjacency_matrix, maxlinewidth)

        plt.text(plotcoordsx[row, split],
                 np.mean([0, plotcoordsy[row, split]]),
                 g1label, fontsize=12)
        
        plt.text(plotcoordsx[row+1, split],
                 np.mean([0, plotcoordsy[row+1, split]]),
                 g2label, fontsize=12)

        plt.text(plotcoordsx[row, split], plotcoordsy[row, split],
                 '     '+str(Z1count), color='red',
                 fontsize=12)
        
        plt.text(plotcoordsx[row+1, split], plotcoordsy[row+1, split],
                 '     '+str(Z2count), color='red',
                 fontsize=12)
        
    else:
        plotcoordsx[row, split] = xrel1 + plotcoordsx[int((row+1)/2), split-1]
        plotcoordsx[row+1, split] = xrel2 + plotcoordsx[int((row+1)/2), split-1]

        plotcoordsy[row, split] = yrel1 + plotcoordsy[int((row+1)/2), split-1]
        plotcoordsy[row+1, split] = yrel2 + plotcoordsy[int((row+1)/2), split-1]

        plot_line(plotcoordsx[int((row+1)/2), split-1],
                  plotcoordsx[row, split],
                  plotcoordsy[int((row+1)/2), split-1],
                  plotcoordsy[row, split],
                  Z1count, adjacency_matrix, maxlinewidth)

        plot_line(plotcoordsx[int((row+1)/2), split-1],
                  plotcoordsx[row+1, split],
                  plotcoordsy[int((row+1)/2), split-1],
                  plotcoordsy[row+1, split],
                  Z2count, adjacency_matrix, maxlinewidth)

        if not np.isnan(plotcoordsx[row, split]):
            plt.text(plotcoordsx[row, split],
                     np.mean([plotcoordsy[int((row+1)/2), split-1],
                     plotcoordsy[row, split]]), g1label, fontsize=12)

            plt.text(plotcoordsx[row, split], plotcoordsy[row, split],
                     '     '+str(Z1count), color='red', fontsize=12)

        if not np.isnan(plotcoordsx[row+1, split]):
            plt.text(plotcoordsx[row+1, split],
                     np.mean([plotcoordsy[int((row+1)/2), split-1],
                     plotcoordsy[row+1, split]]), g2label, fontsize=12)

            plt.text(plotcoordsx[row+1, split], plotcoordsy[row+1, split],
                     '     '+str(Z2count), color='red', fontsize=12)
                    
    return plotcoordsx, plotcoordsy

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage:')
        print('   python3 sCSC.py [adjaency matrix] [num eigvecs] ' +
              '[directory for files]')
        sys.exit(0)
    else:
        adj_mat_file = sys.argv[1]
        numeigvecs = int(sys.argv[2])
        direc = sys.argv[3]

    print('Loading adjacency matrix...')
    try:
        adjacency_matrix = np.load(adj_mat_file)
    except:
        raise RuntimeError('Adjacency matrix file must be a numpy array.')

    # Obtain eigenvectors for clustering and graph Laplacian
    evec, L = get_eigenvectors_for_dendrogram(adjacency_matrix, numeigvecs)

    # Get binary codes for each orthogonal process
    binary = get_binary_codes(evec, numeigvecs)

    # Initialize dendrogram matrices
    groupdist, plotcoordsx, plotcoordsy, index_dict = initialize_dendrogram_matrices(numeigvecs)

    # Loop through eigenvectors, from largest to last one included in analysis
    print('Constructing sCSC dendrogram...')
    for split in np.arange(numeigvecs):
        # Output current dendrogram level
        print("  Split %i of %i" % (split+1, numeigvecs))
        
        # Initialize level of dendrogram
        row = 0
        
        # Loop through binary codes
        for combos in np.arange(0, 2**(split+1), 2):
            g1label, g2label, Z1count, Z2count, Z1, Z2, index_dict = get_subgroup_similarity(
                                                    evec, L, binary, combos, split, index_dict)
            
            # Record dissimilarity metric for group 1
            groupdist[row, split] = Z1
        
            # Record dissimilarity metric for group 2
            groupdist[row+1, split] = Z2
            
            # Plot this level of the dendrogram
            plotcoordsx, plotcoordsy = plot_next_dendrogram_split(
                                       plotcoordsx, plotcoordsy, row, split, adjacency_matrix,
                                       g1label, g2label, Z1count, Z2count, Z1, Z2)
        
            # Advance dendrogram level counter
            row = row + 2

    # Optionally add plot specifications
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    # Make sure directory exists
    try:
        os.stat(direc)
        print('Writing files to existing directory %s...' % direc)
    except:
        os.mkdir(direc)
        print('Creating new directory %s...' % direc)

    # Save figure
    dend_filename = direc  + '/dendrogram.pdf'
    print('Saving dendrogram to %s...' % dend_filename)
    plt.savefig(dend_filename)

    # Export sCSC Dendrogram distance matrix
    dist_filename = direc + '/groupdist.npy'
    print('Saving distance matrix to %s...' % dist_filename)
    np.save(dist_filename, groupdist)

    # Save dictionary of member indices
    index_filename = direc  + '/index_dict.pkl'
    print('Saving index dictionary to %s...' % index_filename)
    with open(index_filename, 'wb') as f:
        pickle.dump(index_dict, f)

    print('Done.')
