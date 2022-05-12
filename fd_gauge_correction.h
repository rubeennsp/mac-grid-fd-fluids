#ifndef FD_GAUGE_CORRECTION_H
#define FD_GAUGE_CORRECTION_H

void fd_node_poisson_solve_3d(
    int ni, int nj, int nk,       // Cell count in each axis (NOT node count)
    double *lap,                  // Laplacian frequencies (ni-1, nj-1, nk-1)
    double *node_freqs_out,       // (ni-1, nj-1, nk-1)
    double cell_sidelength_x = 1,
    double cell_sidelength_y = 1,
    double cell_sidelength_z = 1
);

// Not tested yet
void fd_gauge_correct_3d(
    int ni, int nj, int nk,            // Cell count in each axis (NOT node count)
    double *psi_freqs_x,               // (ni, nj-1, nk-1)
    double *psi_freqs_y,               // (ni-1, nj, nk-1)
    double *psi_freqs_z,               // (ni-1, nj-1, nk)
    double *psi_freqs_x_out = nullptr, // (ni, nj-1, nk-1)
    double *psi_freqs_y_out = nullptr, // (ni-1, nj, nk-1)
    double *psi_freqs_z_out = nullptr, // (ni-1, nj-1, nk)
    double *phi_freqs_out = nullptr,   // (ni-1, nj-1, nk-1)
    double cell_sidelength_x = 1,
    double cell_sidelength_y = 1,          
    double cell_sidelength_z = 1
);

//#############################################################################
// General 3d grid DCT/DST and inverse.
//#############################################################################

enum class Symmetry {
    Even,
    Odd,
};

void to_frequencies_3d(
    const int fluid_grid_size[3], // Cell count in each axis (NOT node count)
    const Symmetry symmetries[3],
    double *data, double *out
);

void from_frequencies_3d(
    const int fluid_grid_size[3], // Cell count in each axis (NOT node count)
    const Symmetry symmetries[3],
    double *data, double *out
);



//#############################################################################
// 3d grid frequency domain transforms and inverses specific to certain fluid data.
//#############################################################################

extern const Symmetry edges_3d_symmetries_x[3];
extern const Symmetry edges_3d_symmetries_y[3];
extern const Symmetry edges_3d_symmetries_z[3];

extern const Symmetry faces_3d_symmetries_x[3];
extern const Symmetry faces_3d_symmetries_y[3];
extern const Symmetry faces_3d_symmetries_z[3];

extern const Symmetry nodes_3d_symmetries[3];

extern const Symmetry cells_3d_symmetries[3];


// Works on interior node values
void to_frequencies_nodes_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *node,                // (ni-1, nj-1, nk-1)
    double *node_out = nullptr   // (ni-1, nj-1, nk-1)
);

// Works on interior node values
void from_frequencies_nodes_3d(
    int ni, int nj, int nk,      // Cell count in each axis (NOT node count)
    double *node,                // (ni-1, nj-1, nk-1)
    double *node_out = nullptr   // (ni-1, nj-1, nk-1)
);


// Works on interior edge values
void to_frequencies_edges_3d(
    int ni, int nj, int nk,        // Cell count in each axis (NOT node count)
    double *edge_x,                // (ni, nj-1, nk-1)
    double *edge_y,                // (ni-1, nj, nk-1)
    double *edge_z,                // (ni-1, nj-1, nk)
    double *edge_x_out = nullptr,  // (ni, nj-1, nk-1)
    double *edge_y_out = nullptr,  // (ni-1, nj, nk-1)
    double *edge_z_out = nullptr   // (ni-1, nj-1, nk)
);

// Works on interior edge values
void from_frequencies_edges_3d(
    int ni, int nj, int nk,        // Cell count in each axis (NOT node count)
    double *edge_x,                // (ni, nj-1, nk-1)
    double *edge_y,                // (ni-1, nj, nk-1)
    double *edge_z,                // (ni-1, nj-1, nk)
    double *edge_x_out = nullptr,  // (ni, nj-1, nk-1)
    double *edge_y_out = nullptr,  // (ni-1, nj, nk-1)
    double *edge_z_out = nullptr   // (ni-1, nj-1, nk)
);


// Works on interior face values
void to_frequencies_faces_3d(
    int ni, int nj, int nk,        // Cell count in each axis (NOT node count)
    double *face_x,                // (ni-1, nj, nk)
    double *face_y,                // (ni, nj-1, nk)
    double *face_z,                // (ni, nj, nk-1)
    double *face_x_out = nullptr,  // (ni-1, nj, nk)
    double *face_y_out = nullptr,  // (ni, nj-1, nk)
    double *face_z_out = nullptr   // (ni, nj, nk-1)
);

// Works on interior face values
void from_frequencies_faces_3d(
    int ni, int nj, int nk,        // Cell count in each axis (NOT node count)
    double *face_x,                // (ni-1, nj, nk)
    double *face_y,                // (ni, nj-1, nk)
    double *face_z,                // (ni, nj, nk-1)
    double *face_x_out = nullptr,  // (ni-1, nj, nk)
    double *face_y_out = nullptr,  // (ni, nj-1, nk)
    double *face_z_out = nullptr   // (ni, nj, nk-1)
);


// Works on interior cell values
void to_frequencies_cells_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *cell,               // (ni, nj, nk)
    double *cell_out            // (ni, nj, nk)
);

// Works on interior cell values
void from_frequencies_cells_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *cell,               // (ni, nj, nk)
    double *cell_out            // (ni, nj, nk)
);

#endif