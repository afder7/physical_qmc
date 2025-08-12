# Physics-inspired low-discrepancy sequence generator for Monte Carlo integration
This branch `second` is about my 2nd semester of working on this research.

Further contents of this `README.md` describe each file which can be found in this directory.


### `generator_with_graph.py`
generates points and builds a graph of discrepancy over time of the simulation. this version of generator chooses the best set of points in terms of discrepancy and and right after the termination of generation compares the discrepancy of the obtained set of points with the set of Sobol points with the same number of points and dimensionality.
