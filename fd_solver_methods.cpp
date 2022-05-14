#include <fftw3.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include "fd_solver_methods.h"


//#############################################################################
// Helpers to deal with odd vs. even signal differences
//#############################################################################


inline int array_size_for_fftw(int num_cells, Symmetry symmetry) {
    // The number of cells along an axis is used as the half-period for the DCT/DST.
    // The illustration below has 4 cells along the x-axis and 2 cells along the y-axis,
    // regardless of whether we're dealing with quantities on the nodes, edges, or faces, and in 2d or 3d.
    //    _ _ _ _
    //   |_|_|_|_| 
    //   |_|_|_|_|
    //
    switch (symmetry) {
    case Symmetry::Even: return num_cells;
    case Symmetry::Odd : return num_cells - 1; // For DST-I, the array length is one lower than the half-period.
    }                                          // This is indeed the correct degree of freedom when considering the 0 boundaries.
    throw 1;
}

/*
    FFTW Documentation for DST/DCT types/kinds
    https://www.fftw.org/fftw3_doc/1d-Real_002dodd-DFTs-_0028DSTs_0029.html
    https://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html
*/

inline fftw_r2r_kind r2r_kind_from_symmetry(Symmetry symmetry) {
    switch (symmetry) {
    case Symmetry::Even: return FFTW_REDFT10;
    case Symmetry::Odd : return FFTW_RODFT00;
    }
    throw 1;
}

inline fftw_r2r_kind inverse_r2r_kind_from_symmetry(Symmetry symmetry) {
    switch (symmetry) {
    case Symmetry::Even: return FFTW_REDFT01;
    case Symmetry::Odd : return FFTW_RODFT00;
    }
    throw 1;
}



//#############################################################################
// General 3d grid DCT/DST and inverse.
//#############################################################################


// Works on interior psi values
void to_frequencies_3d(
    const int fluid_grid_size[3], // Cell count in each axis (NOT node count)
    const Symmetry symmetries[3],
    double *data, double *out
) {
    // Perform in-place if not provided with output location
    if (!out) out = data;

    // Build plan with the right sizes

    const int freq_grid_size[3] = {
        array_size_for_fftw(fluid_grid_size[0], symmetries[0]),
        array_size_for_fftw(fluid_grid_size[1], symmetries[1]),
        array_size_for_fftw(fluid_grid_size[2], symmetries[2]),
    };

    const fftw_r2r_kind transform_kinds[3] = {
        r2r_kind_from_symmetry(symmetries[0]),
        r2r_kind_from_symmetry(symmetries[1]),
        r2r_kind_from_symmetry(symmetries[2]),
    };

    fftw_plan plan = fftw_plan_r2r_3d(
        freq_grid_size[0], freq_grid_size[1], freq_grid_size[2],
        data, out,
        transform_kinds[0], transform_kinds[1], transform_kinds[2],
        FFTW_ESTIMATE
    );

    // Execute
    fftw_execute(plan);

    // Cleanup
    fftw_destroy_plan(plan);
}


void from_frequencies_3d(
    const int fluid_grid_size[3], // Cell count in each axis (NOT node count)
    const Symmetry symmetries[3],
    double *data, double *out
) {
    // Perform in-place if not provided with output location
    if (!out) out = data;

    // Build plan with the right sizes

    const int freq_grid_size[3] = {
        array_size_for_fftw(fluid_grid_size[0], symmetries[0]),
        array_size_for_fftw(fluid_grid_size[1], symmetries[1]),
        array_size_for_fftw(fluid_grid_size[2], symmetries[2]),
    };

    const fftw_r2r_kind transform_kinds[3] = {
        inverse_r2r_kind_from_symmetry(symmetries[0]),
        inverse_r2r_kind_from_symmetry(symmetries[1]),
        inverse_r2r_kind_from_symmetry(symmetries[2]),
    };

    fftw_plan plan = fftw_plan_r2r_3d(
        freq_grid_size[0], freq_grid_size[1], freq_grid_size[2],
        data, out,
        transform_kinds[0], transform_kinds[1], transform_kinds[2],
        FFTW_ESTIMATE
    );

    // Execute
    fftw_execute(plan);

    // Cleanup
    fftw_destroy_plan(plan);

    // IDST/IDCT needs to be downscaled by the full-period (ni, nj, nk are half-periods)
    double multiplier = 1. / double(8 * fluid_grid_size[0] * fluid_grid_size[1] * fluid_grid_size[2]);
    int num_freq_grid_entries = freq_grid_size[0] * freq_grid_size[1] * freq_grid_size[2];
    for (int i = 0; i < num_freq_grid_entries; i++) {
        out[i] *= multiplier;
    }
}



//#############################################################################
// Solvers in frequency domain
//#############################################################################

// Private vector helpers
namespace {
    struct vec3 { double x, y, z; };

    inline double dot(vec3 v, vec3 w) {
        return v.x * w.x + v.y * w.y + v.z * w.z;
    }

    inline double mag2(vec3 v) {
        return dot(v, v);
    }

    inline vec3 cross(vec3 v, vec3 w) {
        return {
            v.y * w.z - v.z * w.y,
            v.z * w.x - v.x * w.z,
            v.x * w.y - v.y * w.x,
        };
    }

    inline vec3 operator* (double c, vec3 v) {
        return {
            c * v.x,
            c * v.y,
            c * v.z,
        };
    }

    inline vec3 operator* (vec3 v, double c) {
        return c * v;
    }

    inline vec3 operator+ (vec3 v, vec3 w) {
        return {
            v.x + w.x,
            v.y + w.y,
            v.z + w.z,
        };
    }

    inline vec3 operator- (vec3 v, vec3 w) {
        return {
            v.x - w.x,
            v.y - w.y,
            v.z - w.z,
        };
    }

    inline vec3 operator- (vec3 v) {
        return { -v.x, -v.y, -v.z, };
    }
}


// Helper for indexing flat grids
static const int INVALID_IDX = -1;
inline int grid_index_3d(int ix, int iy, int iz, int lenx, int leny, int lenz) {
    if (
        ix >= 0 || iy >= 0 || iz >= 0
        || ix < lenx || iy < leny || iz < lenz
    )
        return (ix * leny * lenz) + (iy * lenz) + iz;

    return INVALID_IDX;
}


// TODO: Explain what these diff multipliers are
std::vector<double> fd_finite_differencing_multipliers_s2c(int n, double length_between_samples = 1) {
    std::vector<double> v;
    v.reserve(n);

    for (int i = 0; i < n; i++) 
        v.push_back(
            2 * sin(double(i) * M_PI_2 / double(n))
            / length_between_samples
        );

    return v;
}


void fd_node_poisson_solve_3d(
    int ni, int nj, int nk,       // Cell count in each axis (NOT node count)
    double *lap_freqs,            // Laplacian frequencies (ni-1, nj-1, nk-1)
    double *node_freqs_out,       // (ni-1, nj-1, nk-1)
    double cell_sidelength_x,
    double cell_sidelength_y,
    double cell_sidelength_z
) {
    // Indexing
    auto get_idx = [=](int i, int j, int k) { return grid_index_3d(i-1, j-1, k-1, ni-1, nj-1, nk-1); };

    // Modify input array in-place if no output location is provided.
    if (!node_freqs_out) node_freqs_out = lap_freqs;

    // WARNING: If you modify the following, it has to take into account
    // that the output location `node_freqs_out` might be the same as
    // the input location `lap_freqs`.
    // That is, do not read after writing.

    // Finite differencing multipliers
    std::vector<double> diff_s2c_x = fd_finite_differencing_multipliers_s2c(ni, cell_sidelength_x);
    std::vector<double> diff_s2c_y = fd_finite_differencing_multipliers_s2c(nj, cell_sidelength_y);
    std::vector<double> diff_s2c_z = fd_finite_differencing_multipliers_s2c(nk, cell_sidelength_z);

    // The for loops here can start at 1 because all axes use DST-I and don't have the 0 frequency.
    for (int i = 1; i < ni; i++) {
        for (int j = 1; j < nj; j++) {
            for (int k = 1; k < nk; k++) {
                int node_idx = get_idx(i, j, k);
                assert(node_idx != INVALID_IDX); // True because the for loop ranges are correct. Assert anyways.

                // TODO: Explain what these diff vectors are
                vec3 diff_s2c_ijk {
                    diff_s2c_x[i],
                    diff_s2c_y[j],
                    diff_s2c_z[k],
                };
                vec3 diff_c2s_ijk = -diff_s2c_ijk;

                // In this frequency domain, taking the laplacian of f is just dot(diff, diff) * f
                // although we just need to mind which diff multipliers are used. (Sine-to-cosine or the other way around.)
                // If g is the laplacian of f, then f = g / dot(diff, diff)
                node_freqs_out[node_idx] = lap_freqs[node_idx] / dot(diff_c2s_ijk, diff_s2c_ijk);
                                                               // equiv = -mag2(diff)
            }
        }
    }
}


void fd_cell_poisson_solve_3d(
    int ni, int nj, int nk,       // Cell count in each axis (NOT node count)
    double *lap_freqs,            // Laplacian frequencies (ni, nj, nk)
    double *cell_freqs_out,       // (ni, nj, nk)
    double cell_sidelength_x,
    double cell_sidelength_y,
    double cell_sidelength_z
) {
    // Indexing
    auto get_idx = [=](int i, int j, int k) { return grid_index_3d(i, j, k, ni, nj, nk); };

    // Modify input array in-place if no output location is provided.
    if (!cell_freqs_out) cell_freqs_out = lap_freqs;

    // WARNING: If you modify the following, it has to take into account
    // that the output location `cell_freqs_out` might be the same as
    // the input location `lap_freqs`.
    // That is, do not read after writing.

    // Finite differencing multipliers
    std::vector<double> diff_s2c_x = fd_finite_differencing_multipliers_s2c(ni, cell_sidelength_x);
    std::vector<double> diff_s2c_y = fd_finite_differencing_multipliers_s2c(nj, cell_sidelength_y);
    std::vector<double> diff_s2c_z = fd_finite_differencing_multipliers_s2c(nk, cell_sidelength_z);

    // The for loops here can start at 1 because all axes use DST-I and don't have the 0 frequency.
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                // Special case: The "constant shift" frequency mode can take on any value
                if (i == 0 && j == 0 && k == 0) {
                    cell_freqs_out[0] = 0; // set to zero, although any other value is OK as a constant shift that doesn't change the laplacian
                    continue;
                }

                int cell_idx = get_idx(i, j, k);
                assert(cell_idx != INVALID_IDX); // True because the for loop ranges are correct. Assert anyways.

                // TODO: Explain what these diff vectors are
                vec3 diff_s2c_ijk {
                    diff_s2c_x[i],
                    diff_s2c_y[j],
                    diff_s2c_z[k],
                };
                vec3 diff_c2s_ijk = -diff_s2c_ijk;

                // In this frequency domain, taking the laplacian of f is just dot(diff, diff) * f
                // although we just need to mind which diff multipliers are used. (Sine-to-cosine or the other way around.)
                // If g is the laplacian of f, then f = g / dot(diff, diff)
                cell_freqs_out[cell_idx] = lap_freqs[cell_idx] / dot(diff_c2s_ijk, diff_s2c_ijk);
                                                               // equiv = -mag2(diff)
            }
        }
    }
}


// Not tested yet. Probably better to not use this because
// it requires more data to have gone through frequency transforms
// (three grids instead of just the one needed by fd_node_poisson_solve_3d)
void fd_gauge_correct_3d(
    int ni, int nj, int nk,            // Cell count in each axis (NOT node count)
    double *psi_freqs_x,               // (ni, nj-1, nk-1)
    double *psi_freqs_y,               // (ni-1, nj, nk-1)
    double *psi_freqs_z,               // (ni-1, nj-1, nk)
    double *psi_freqs_x_out,           // (ni, nj-1, nk-1)
    double *psi_freqs_y_out,           // (ni-1, nj, nk-1)
    double *psi_freqs_z_out,           // (ni-1, nj-1, nk)
    double *phi_freqs_out,             // (ni-1, nj-1, nk-1)
    double cell_sidelength_x,
    double cell_sidelength_y,          
    double cell_sidelength_z
) {
    // Indexing
    auto psi_x_idx = [=](int i, int j, int k) { return grid_index_3d(i  , j-1, k-1, ni  , nj-1, nk-1); };
    auto psi_y_idx = [=](int i, int j, int k) { return grid_index_3d(i-1, j  , k-1, ni-1, nj  , nk-1); };
    auto psi_z_idx = [=](int i, int j, int k) { return grid_index_3d(i-1, j-1, k  , ni-1, nj-1, nk  ); };

    auto phi_idx   = [=](int i, int j, int k) { return grid_index_3d(i-1, j-1, k-1, ni-1, nj-1, nk-1); };

    // Modify input array in-place if no output location is provided.
    if (!psi_freqs_x_out) psi_freqs_x_out = psi_freqs_x;
    if (!psi_freqs_y_out) psi_freqs_y_out = psi_freqs_y;
    if (!psi_freqs_z_out) psi_freqs_z_out = psi_freqs_z;

    // Finite differencing multipliers
    std::vector<double> diff_s2c_x = fd_finite_differencing_multipliers_s2c(ni, cell_sidelength_x);
    std::vector<double> diff_s2c_y = fd_finite_differencing_multipliers_s2c(nj, cell_sidelength_y);
    std::vector<double> diff_s2c_z = fd_finite_differencing_multipliers_s2c(nk, cell_sidelength_z);

    // WARNING: If you modify the following, it has to take into account
    // that the output location `psi_freqs_x_out` might be the same as
    // the input location `psi_freqs_x`.
    // That is, do not read after writing.

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                int psi_x_idx_ijk = psi_x_idx(i, j, k);
                int psi_y_idx_ijk = psi_y_idx(i, j, k);
                int psi_z_idx_ijk = psi_z_idx(i, j, k);

                vec3 psi_ijk {
                    psi_x_idx_ijk == INVALID_IDX ? 0 : psi_freqs_x[psi_x_idx_ijk],
                    psi_y_idx_ijk == INVALID_IDX ? 0 : psi_freqs_y[psi_y_idx_ijk],
                    psi_z_idx_ijk == INVALID_IDX ? 0 : psi_freqs_z[psi_z_idx_ijk],
                };

                vec3 diff_s2c_ijk {
                    diff_s2c_x[i],
                    diff_s2c_y[j],
                    diff_s2c_z[k],
                };
                vec3 diff_c2s_ijk = - diff_s2c_ijk;

                // Laplacian of phi = phi * dot(diff_s2c, diff_c2s) or equivalently phi * (-mag(diff))
                // Divergence of psi = dot(diff, psi)
                // Equating the two lines above gives the following solution for phi
                double phi_ijk = dot(diff_c2s_ijk, psi_ijk) / dot(diff_c2s_ijk, diff_s2c_ijk);

                // In the frequency domain, the gradient of phi is just phi * diff. (scalar * vector)
                vec3 grad_phi_ijk = phi_ijk * diff_s2c_ijk;

                // Divergence-free psi
                vec3 corrected_psi_ijk = psi_ijk - grad_phi_ijk;

                // Write out to psi
                if (psi_x_idx_ijk != INVALID_IDX) psi_freqs_x_out[psi_x_idx_ijk] = corrected_psi_ijk.x;
                if (psi_y_idx_ijk != INVALID_IDX) psi_freqs_y_out[psi_y_idx_ijk] = corrected_psi_ijk.y;
                if (psi_z_idx_ijk != INVALID_IDX) psi_freqs_z_out[psi_z_idx_ijk] = corrected_psi_ijk.z;
                
                if (phi_freqs_out) {
                    int phi_idx_ijk = phi_idx(i, j, k);
                    if (phi_idx_ijk != INVALID_IDX) phi_freqs_out[phi_idx_ijk] = phi_ijk;
                }
            }
        }
    }
}



//#############################################################################
// 3d grid frequency domain transforms and inverses specific to certain elements of the staggered grid.
// 
// Maybe we don't need every single one of these but they're a good reference
// to help understand what's going on with the odd/even symmetries anyways.
//#############################################################################

const Symmetry edges_3d_symmetries_x[3] = { Symmetry::Even, Symmetry::Odd, Symmetry::Odd, };
const Symmetry edges_3d_symmetries_y[3] = { Symmetry::Odd, Symmetry::Even, Symmetry::Odd, };
const Symmetry edges_3d_symmetries_z[3] = { Symmetry::Odd, Symmetry::Odd, Symmetry::Even, };

const Symmetry faces_3d_symmetries_x[3] = { Symmetry::Odd, Symmetry::Even, Symmetry::Even, };
const Symmetry faces_3d_symmetries_y[3] = { Symmetry::Even, Symmetry::Odd, Symmetry::Even, };
const Symmetry faces_3d_symmetries_z[3] = { Symmetry::Even, Symmetry::Even, Symmetry::Odd, };

const Symmetry nodes_3d_symmetries[3] = { Symmetry::Odd, Symmetry::Odd, Symmetry::Odd };

const Symmetry cells_3d_symmetries[3] = { Symmetry::Even, Symmetry::Even, Symmetry::Even, };


// Works on interior edge values
void to_frequencies_edges_3d(
    int ni, int nj, int nk,
    double *edge_x,                // (ni, nj-1, nk-1)
    double *edge_y,                // (ni-1, nj, nk-1)
    double *edge_z,                // (ni-1, nj-1, nk)
    double *edge_x_out,            // (ni, nj-1, nk-1)
    double *edge_y_out,            // (ni-1, nj, nk-1)
    double *edge_z_out             // (ni-1, nj-1, nk)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    to_frequencies_3d(fluid_grid_size, edges_3d_symmetries_x, edge_x, edge_x_out);
    to_frequencies_3d(fluid_grid_size, edges_3d_symmetries_y, edge_y, edge_y_out);
    to_frequencies_3d(fluid_grid_size, edges_3d_symmetries_z, edge_z, edge_z_out);
}


// Works on interior edge values
void from_frequencies_edges_3d(
    int ni, int nj, int nk,
    double *edge_x,                // (ni, nj-1, nk-1)
    double *edge_y,                // (ni-1, nj, nk-1)
    double *edge_z,                // (ni-1, nj-1, nk)
    double *edge_x_out,            // (ni, nj-1, nk-1)
    double *edge_y_out,            // (ni-1, nj, nk-1)
    double *edge_z_out             // (ni-1, nj-1, nk)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    from_frequencies_3d(fluid_grid_size, edges_3d_symmetries_x, edge_x, edge_x_out);
    from_frequencies_3d(fluid_grid_size, edges_3d_symmetries_y, edge_y, edge_y_out);
    from_frequencies_3d(fluid_grid_size, edges_3d_symmetries_z, edge_z, edge_z_out);
}


// Works on interior face values
void to_frequencies_faces_3d(
    int ni, int nj, int nk,
    double *face_x,                // (ni-1, nj, nk)
    double *face_y,                // (ni, nj-1, nk)
    double *face_z,                // (ni, nj, nk-1)
    double *face_x_out,            // (ni-1, nj, nk)
    double *face_y_out,            // (ni, nj-1, nk)
    double *face_z_out             // (ni, nj, nk-1)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    to_frequencies_3d(fluid_grid_size, faces_3d_symmetries_x, face_x, face_x_out);
    to_frequencies_3d(fluid_grid_size, faces_3d_symmetries_y, face_y, face_y_out);
    to_frequencies_3d(fluid_grid_size, faces_3d_symmetries_z, face_z, face_z_out);
}


// Works on interior face values
void from_frequencies_faces_3d(
    int ni, int nj, int nk,
    double *face_x,                // (ni-1, nj, nk)
    double *face_y,                // (ni, nj-1, nk)
    double *face_z,                // (ni, nj, nk-1)
    double *face_x_out,            // (ni-1, nj, nk)
    double *face_y_out,            // (ni, nj-1, nk)
    double *face_z_out             // (ni, nj, nk-1)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    from_frequencies_3d(fluid_grid_size, faces_3d_symmetries_x, face_x, face_x_out);
    from_frequencies_3d(fluid_grid_size, faces_3d_symmetries_y, face_y, face_y_out);
    from_frequencies_3d(fluid_grid_size, faces_3d_symmetries_z, face_z, face_z_out);
}


// Works on interior node values
void to_frequencies_nodes_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *node,               // (ni-1, nj-1, nk-1)
    double *node_out            // (ni-1, nj-1, nk-1)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    to_frequencies_3d(fluid_grid_size, nodes_3d_symmetries, node, node_out);
}


// Works on interior node values
void from_frequencies_nodes_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *node,                // (ni-1, nj-1, nk-1)
    double *node_out             // (ni-1, nj-1, nk-1)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    from_frequencies_3d(fluid_grid_size, nodes_3d_symmetries, node, node_out);
}


// Works on interior cell values
void to_frequencies_cells_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *cell,               // (ni, nj, nk)
    double *cell_out = nullptr  // (ni, nj, nk)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    to_frequencies_3d(fluid_grid_size, cells_3d_symmetries, cell, cell_out);
}


// Works on interior cell values
void from_frequencies_cells_3d(
    int ni, int nj, int nk,     // Cell count in each axis (NOT node count)
    double *cell,                // (ni, nj, nk)
    double *cell_out = nullptr   // (ni, nj, nk)
) {
    int fluid_grid_size[3] = {ni, nj, nk};
    from_frequencies_3d(fluid_grid_size, cells_3d_symmetries, cell, cell_out);
}