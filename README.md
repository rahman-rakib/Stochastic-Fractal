# Stochastic-Fractal

# Folders
## data
`data` folder contains all the data generated from simulation.

## figures
all the figures generated from data in `data` folder or from analytic functions

## main_jupyter
contains jupyter notebook, i.e. `.ipynb`, files that uses class structured code for
 simulation and analytic solutions


## main_py
pure python codes written using object oriented style, i.e. python class.
### `stochastic_fragmentation*.py` structure is the simulation code

Sl | class name | super class | description
---|------------|-------------|----------------
1. | StochasticFragmentation | - | super class for all other classes in this folder. Basic iteration step is defined here.
2. | NumberLength | StochasticFragmentation | to generate ensemble average of 0-th and 1st moment.
3. | Moment       | StochasticFragmentation | finding ensemble average of the n-th moment.
4. | TrueLengths  | StochasticFragmentation | generate a list of segment lengths for a number of realization. for data collapse plot.

### `analytic_solution.py` is for analytic solution
for $\alpha=1,2,3$ the analytic solution is known. 
and fractal dimension $d_f$ can be found for any $\alpha$ and $p$ 

## notebooks
all codes written in jupyter notebook.
can generate all necessary data and plots.

## `.sh` files
we have used it to open jupyter notebook in this folder, where `.sh` files are.


## `stochastic_fractal_m.py` or the Fragmentation module
This module must be run whenever one tries to run the classes in jupyter notebook or external python file.
