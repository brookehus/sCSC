Data
====

Requirements:
+ numpy

This directory contains the three fluid mechanics datasets used in this work (`quad_gyre`, `bickley_jet`, `vortex`) and the and the four adjacency matrices analyzed in this work (all the previous and `protein_g`).

`*_trajectories.npy` files have three dimensions: (0) the number of dimensions in the dataset, (1) the number of trajectories in the dataset, (2) the number of time points in each trajectory.

`*_adjmat.npy` files are symmetric square matrices where the dimension is the number of trajectories and each entry is the dissimilarity measure between the row and column indices.

Numpy files can be loaded in python using:

```my_data = np.load('my_filename.npy')```

The dimensions are shown by running:

```my_data.shape```

Note
----

Please note that the vortex data trajectectory 416 (0-indexed) has NaNs. The corresponding adjacency matrix is not affected.
