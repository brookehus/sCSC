Example python script and notebooks
===================================

Requirements:
+ numpy
+ matplotlib
+ scipy
+ pickle (optional; for saving dictionaries)

For notebooks:
+ Jupyter

The calculation requires a square adjacency matrix in numpy format.

The python script `sCSC.py` can be run from the command line using:

```python3 sCSC.py [adjaency matrix] [num eigvecs] [directory for files]```

The `sCSC_dendrogram_creation` notebook performs the same calculation as the `sCSC.py` in an interactive Jupyter environment.

The `sCSC_data_overview` notebook briefly overviews the format of the data contained in the `sCSC/data` directory.
