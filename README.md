# Physics-inspired low-discrepancy sequence generator for Monte Carlo integration
This branch `first` and `README.md` file specifically represents all files created and used during my 1st semester (2nd in general at MIPT (Moscow Institute of Physics and Technology)) of work on this research.

Report on this research could be found in `4.pdf` file in this branch.

My current work and new code for the intended library for QMC Integration could be found in the branch `second`.

Further contents of this `README.md` describe each file which can be found in this directory.


### `discrepancy.py`
change num, dim to get l2 of each file of these num and dim

### `electric_model.py`
2d generator without discrepancy decrease graph output

### `error_calc_generator.py`
build graph with error of calculation depending on the number of points
used in integral calculation

### `error_graph.py`
bad no use

### `generator.py`
v3 point generator

### `gravitational_model.py`
old point generator

### `integrals.py`
tests of multiple integrals error of sobol and our points

### `projections.py`
builds graph of projections of generated point sequences on 2d plane to
check if the projections are uniform (should be rather checked as a
graph over iterations which represents discrepancy change on several
planes