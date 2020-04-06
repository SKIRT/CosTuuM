Scripts used to generate the figures in the CosTuuM code paper.

The scripts have been split into two categories: scripts that run 
CosTuuM simulations and require a lot of computational resources, and 
scripts that create figures. Data from the former is communicated to the 
latter through intermediary binary files written using numpy.tofile and 
read using numpy.fromfile.

The dependency graph for the scripts is contained in the Makeflow 
workflow file `create_figures.makeflow`.
