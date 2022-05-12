Compilation in Ruben's Linux machine: 
```
g++ -c fd_solver_methods.cpp
g++ main.cpp fd_solver_methods.o -lfftw3 -lm -o main
```

@JC, the poisson solve on interior grid nodes is `fd_node_poisson_solve_3d`, declared in `fd_solver_methods.h`. That file also contains the methods to convert between spatial-domain data and frequency-domain data. I didn't take the time to figure out how to customize stride/padding/row-order in FFTW, so these methods only accept flat compact arrays of the right sizes that FFTW would typically receive. They want C-order (row-major), we want to pass in _only_ the interior node grid, so `grid(i, j, k) = grid[i*(nj-1)*(nk-1) + j*(nk-1) + k]`. You might have to copy over values from your grid to fit the size and order requirements.

Conversion methods like `to_frequencies_nodes_3d` have the expected input array sizes annotated relative to the fluid grid size. As a convention, the methods here always expect the `(ni, nj, nk)` to be specified as the cell-count of the fluid grid and never the node-count or edge-count or anything else.

The poisson solve assumes zero at the boundary nodes for both the output and the laplacian (with appropriate ghost values outside the rectangular boundary).

Please check out the example usage of `fd_node_poisson_solve_3d`, `to_frequencies_nodes_3d`, and `from_frequencies_nodes_3d` within the `test2` function in `main.cpp`.
I hope this is enough, let me know if you have any questions!
