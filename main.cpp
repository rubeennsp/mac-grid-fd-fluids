#include "fd_solver_methods.h"
#include <iostream>
#include <string>
#include <chrono>

void print_title(std::string title) {
    std::cout
        << std::string(80, '#') << std::endl
        << "# " << title << std::endl
        << std::string(80, '#') << std::endl
        << std::endl;
}

void print_grid_3d(int ni, int nj, int nk, double *grid) {
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                std::cout << grid[i*nj*nk + j*nk + k] << " ";
            }
            std::cout << std::endl;
        }
        std:: cout << std::endl;
    }
}

void print_grid_2d(int ni, int nj, double *grid) {
    return print_grid_3d(1, ni, nj, grid);
}

void test1() {
    print_title("Test 1: Frequency transform on nodes");

    double grid[3][5][7] = {
        {
            {111, 112, 113, 114, 115, 116, 117,},
            {121, 122, 123, 124, 125, 126, 127,},
            {131, 132, 133, 134, 135, 136, 137,},
            {141, 142, 143, 144, 145, 146, 147,},
            {151, 152, 153, 154, 155, 156, 157,},
        },
        {
            {211, 212, 213, 214, 215, 216, 217,},
            {221, 222, 223, 224, 225, 226, 227,},
            {231, 232, 233, 234, 235, 236, 237,},
            {241, 242, 243, 244, 245, 246, 247,},
            {251, 252, 253, 254, 255, 256, 257,},
        },
        {
            {311, 312, 313, 314, 315, 316, 317,},
            {321, 322, 323, 324, 325, 326, 327,},
            {331, 332, 333, 334, 335, 336, 337,},
            {341, 342, 343, 344, 345, 346, 347,},
            {351, 352, 353, 354, 355, 356, 357,},
        },
    };

    std::cout << "Original" << std::endl;
    print_grid_3d(3, 5, 7, &grid[0][0][0]);
    std::cout << std::endl;

    // Convert into frequencies in place
    to_frequencies_nodes_3d(4, 6, 8, &grid[0][0][0]);

    std::cout << "Frequencies" << std::endl;
    print_grid_3d(3, 5, 7, &grid[0][0][0]);
    std::cout << std::endl;

    // Convert back into original data in place
    from_frequencies_nodes_3d(4, 6, 8, &grid[0][0][0]);

    std::cout << "Reconstructed" << std::endl;
    print_grid_3d(3, 5, 7, &grid[0][0][0]);
    std::cout << std::endl;
}

void test2() {
    print_title("Test 2: Frequency-domain poisson solve");

    // In a (4, 6, 8) fluid grid, the interior nodes are a (3, 5, 7) grid.
    double target_laplacian[3][5][7] =
    {
        {
            {111, 112, 113, 114, 115, 116, 117,},
            {121, 122, 123, 124, 125, 126, 127,},
            {131, 132, 133, 134, 135, 136, 137,},
            {141, 142, 143, 144, 145, 146, 147,},
            {151, 152, 153, 154, 155, 156, 157,},
        },
        {
            {211, 212, 213, 214, 215, 216, 217,},
            {221, 222, 223, 224, 225, 226, 227,},
            {231, 232, 233, 234, 235, 236, 237,},
            {241, 242, 243, 244, 245, 246, 247,},
            {251, 252, 253, 254, 255, 256, 257,},
        },
        {
            {311, 312, 313, 314, 315, 316, 317,},
            {321, 322, 323, 324, 325, 326, 327,},
            {331, 332, 333, 334, 335, 336, 337,},
            {341, 342, 343, 344, 345, 346, 347,},
            {351, 352, 353, 354, 355, 356, 357,},
        },
    };


    double sidelengths[3] = { 1, 2, 4 }; // Can try other sidelengths

    std::cout << "Target laplacian (at interior nodes. Boundary nodes are assumed to be 0)" << std::endl;
    print_grid_3d(3, 5, 7, &target_laplacian[0][0][0]);
    std::cout << std::endl;

    // Convert into frequencies in-place.
    // In a (4, 6, 8) fluid grid, the interior nodes are a (3, 5, 7) grid.
    to_frequencies_nodes_3d(4, 6, 8, &target_laplacian[0][0][0]);

    // Poisson solve in frequency domain
    double result[3][5][7];
    fd_node_poisson_solve_3d(4, 6, 8, &target_laplacian[0][0][0], &result[0][0][0], sidelengths[0], sidelengths[1], sidelengths[2]);

    // Convert results back from frequencies to actual node values (in-place)
    from_frequencies_nodes_3d(4, 6, 8, &result[0][0][0]);

    std::cout << "Result" << std::endl;
    print_grid_3d(3, 5, 7, &result[0][0][0]);
    std::cout << std::endl;

    // Calculate the laplacian of the solver output to print for comparison with the target laplacian
    double result_laplacian[3][5][7];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 7; k++) {
                double left   = (i > 0   ? result[i-1][j][k] : 0);
                double right  = (i < 3-1 ? result[i+1][j][k] : 0);
                double bottom = (j > 0   ? result[i][j-1][k] : 0);
                double top    = (j < 5-1 ? result[i][j+1][k] : 0);
                double front  = (k > 0   ? result[i][j][k-1] : 0);
                double back   = (k < 7-1 ? result[i][j][k+1] : 0);
                double center = result[i][j][k];

                result_laplacian[i][j][k] = (
                    ((right + left - 2 * center)) / (sidelengths[0] * sidelengths[0])
                  + ((top + bottom - 2 * center)) / (sidelengths[1] * sidelengths[1])
                  + ((front + back - 2 * center)) / (sidelengths[2] * sidelengths[2])
                );
            }
        }
    }

    std::cout << "Result laplacian" << std::endl;
    print_grid_3d(3, 5, 7, &result_laplacian[0][0][0]);
    std::cout << std::endl;

    std::cout
        << "Please verify that the result laplacian matches the target laplacian." << std::endl
        << "This laplacian was calculated using sidelengths "
            << sidelengths[0] << ", " << sidelengths[1] << ", and " << sidelengths[2] << " "
            << "for the x, y, and z axes respectively." << std::endl
        << "Please refer to the test code to verify that the laplacian is calculated as expected." << std::endl
        << std::endl
        << std::endl;
}

void test3() {
    print_title("Test 3: Basic benchmark for poisson solve");

    int ni = 128, nj = 256, nk = 512;

    double *target_laplacian = new double[(ni-1)*(nj-1)*(nk-1)];

    double sidelengths[3] = { 1, 2, 4 };

    int num_repetitions = 10;
    int print_every = 1;

    std::cout
        << "Provided input:" << std::endl
        << "  Target laplacian at interior nodes, already arranged in a compact array." << std::endl
        << "    (Conceptual) Fluid grid size = (" << ni << ", " << nj << ", " << nk << ")" << std::endl
        << "    Interior node grid size = (" << ni-1 << ", " << nj-1 << ", " << nk-1 << ")" << std::endl
        << "Memory:" << std::endl
        << "  Some pre-allocated grids for (intermediate and final) results." << std::endl
        << "    See code to check where in-place operations happen." << std::endl
        << "Task (to be repeated " << num_repetitions << " times):" << std::endl
        << "  1. Transform target laplacian to frequencies." << std::endl
        << "  2. Perform poisson solve in frequency domain." << std::endl
        << "  3. Convert solve result back from frequency domain into spatial." << std::endl
        << std::endl;

    double *result = new double[(ni-1)*(nj-1)*(nk-1)];
    // Can allocate more arrays and use them to see how non-in-place operations perform.
    // double *result2 = new double[(ni-1)*(nj-1)*(nk-1)];
    // double *result3 = new double[(ni-1)*(nj-1)*(nk-1)];

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_repetitions; i++) {
        if (i % print_every == 0 && i > 0) std::cout << i << " repetitions done" << std::endl << std::flush;

        to_frequencies_nodes_3d(ni, nj, nk, target_laplacian, result); // Write out to result array
        fd_node_poisson_solve_3d(ni, nj, nk, result, result, sidelengths[0], sidelengths[1], sidelengths[2]); // in-place on the result array
        from_frequencies_nodes_3d(ni, nj, nk, result, result); // in-place on the result array
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << std::endl << "All " << num_repetitions << " repetitions done!" << std::endl;
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout
        << "Total elapsed time: " << float_ms.count() / 1000 << " seconds." << std::endl
        << "That is " << float_ms.count() / 1000 / double(num_repetitions) << " seconds per repetition." << std::endl
        << std::endl
        << std::endl;

    delete[] target_laplacian;
    delete[] result;
}

void test4() {
    print_title("Test 4: Frequency domain poisson solve with non-zero dirichlet boundary conditions");

    // In a (4, 6, 8) fluid grid, the interior nodes are a (3, 5, 7) grid.
    double target_laplacian[3][5][7] =
    {
        {
            {111, 112, 113, 114, 115, 116, 117,},
            {121, 122, 123, 124, 125, 126, 127,},
            {131, 132, 133, 134, 135, 136, 137,},
            {141, 142, 143, 144, 145, 146, 147,},
            {151, 152, 153, 154, 155, 156, 157,},
        },
        {
            {211, 212, 213, 214, 215, 216, 217,},
            {221, 222, 223, 224, 225, 226, 227,},
            {231, 232, 233, 234, 235, 236, 237,},
            {241, 242, 243, 244, 245, 246, 247,},
            {251, 252, 253, 254, 255, 256, 257,},
        },
        {
            {311, 312, 313, 314, 315, 316, 317,},
            {321, 322, 323, 324, 325, 326, 327,},
            {331, 332, 333, 334, 335, 336, 337,},
            {341, 342, 343, 344, 345, 346, 347,},
            {351, 352, 353, 354, 355, 356, 357,},
        },
    };

    // Define Dirichlet boundary conditions

    double dirichlet_x_lo[5][7] =
    {
        {411, 412, 413, 414, 415, 416, 417,},
        {421, 422, 423, 424, 425, 426, 427,},
        {431, 432, 433, 434, 435, 436, 437,},
        {441, 442, 443, 444, 445, 446, 447,},
        {451, 452, 453, 454, 455, 456, 457,},
    };
    double dirichlet_x_hi[5][7] =
    {
        {711, 712, 713, 714, 715, 716, 717,},
        {721, 722, 723, 724, 725, 726, 727,},
        {731, 732, 733, 734, 735, 736, 737,},
        {741, 742, 743, 744, 745, 746, 747,},
        {751, 752, 753, 754, 755, 756, 757,},
    };
    double dirichlet_y_lo[3][7] =
    {
        {511, 512, 513, 514, 515, 516, 517,},
        {521, 522, 523, 524, 525, 526, 527,},
        {531, 532, 533, 534, 535, 536, 537,},
    };
    double dirichlet_y_hi[3][7] =
    {
        {811, 812, 813, 814, 815, 816, 817,},
        {821, 822, 823, 824, 825, 826, 827,},
        {831, 832, 833, 834, 835, 836, 837,},
    };
    double dirichlet_z_lo[3][5] =
    {
        {611, 612, 613, 614, 615,},
        {621, 622, 623, 624, 625,},
        {631, 632, 633, 634, 635,},
    };
    double dirichlet_z_hi[3][5] =
    {
        {911, 912, 913, 914, 915,},
        {921, 922, 923, 924, 925,},
        {931, 932, 933, 934, 935,},
    };


    double sidelengths[3] = { 1, 2, 4 }; // Can try other sidelengths

    // Print problem data

    std::cout << "Target laplacian (at interior nodes)" << std::endl;
    print_grid_3d(3, 5, 7, &target_laplacian[0][0][0]);
    std::cout << std::endl;

    std::cout << "Dirichlet boundary condition at low x wall" << std::endl;
    print_grid_2d(5, 7, &dirichlet_x_lo[0][0]);

    std::cout << "Dirichlet boundary condition at high x wall" << std::endl;
    print_grid_2d(5, 7, &dirichlet_x_hi[0][0]);

    std::cout << "Dirichlet boundary condition at low y wall" << std::endl;
    print_grid_2d(3, 7, &dirichlet_y_lo[0][0]);

    std::cout << "Dirichlet boundary condition at high y wall" << std::endl;
    print_grid_2d(3, 7, &dirichlet_y_hi[0][0]);

    std::cout << "Dirichlet boundary condition at low z wall" << std::endl;
    print_grid_2d(3, 5, &dirichlet_z_lo[0][0]);

    std::cout << "Dirichlet boundary condition at high z wall" << std::endl;
    print_grid_2d(3, 5, &dirichlet_z_hi[0][0]);

    std::cout << std::endl;

    // Transform problem to have zero Dirichlet boundary conditions
    // by subtracting boundary_val / sidelength^2 from RHS at the
    // interior nodes that neighbor the boundary.

    for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 7; k++) {
            target_laplacian[  0][j][k] -= dirichlet_x_lo[j][k] / (sidelengths[0] * sidelengths[0]);
            target_laplacian[3-1][j][k] -= dirichlet_x_hi[j][k] / (sidelengths[0] * sidelengths[0]);
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 7; k++) {
            target_laplacian[i][  0][k] -= dirichlet_y_lo[i][k] / (sidelengths[1] * sidelengths[1]);
            target_laplacian[i][5-1][k] -= dirichlet_y_hi[i][k] / (sidelengths[1] * sidelengths[1]);
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 5; j++) {
            target_laplacian[i][j][  0] -= dirichlet_z_lo[i][j] / (sidelengths[2] * sidelengths[2]);
            target_laplacian[i][j][7-1] -= dirichlet_z_hi[i][j] / (sidelengths[2] * sidelengths[2]);
        }
    }

    // Continue with the zero-Dirichlet Poisson solve seen in previous tests

    // Convert into frequencies in-place.
    // In a (4, 6, 8) fluid grid, the interior nodes are a (3, 5, 7) grid.
    to_frequencies_nodes_3d(4, 6, 8, &target_laplacian[0][0][0]);

    // Poisson solve in frequency domain
    double result[3][5][7];
    fd_node_poisson_solve_3d(4, 6, 8, &target_laplacian[0][0][0], &result[0][0][0], sidelengths[0], sidelengths[1], sidelengths[2]);

    // Convert results back from frequencies to actual node values (in-place)
    from_frequencies_nodes_3d(4, 6, 8, &result[0][0][0]);

    std::cout << "Result" << std::endl;
    print_grid_3d(3, 5, 7, &result[0][0][0]);
    std::cout << std::endl;

    // Calculate the laplacian of the solver output to print for comparison with the target laplacian
    double result_laplacian[3][5][7];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 7; k++) {
                double left   = (i > 0   ? result[i-1][j][k] : dirichlet_x_lo[j][k]);
                double right  = (i < 3-1 ? result[i+1][j][k] : dirichlet_x_hi[j][k]);
                double bottom = (j > 0   ? result[i][j-1][k] : dirichlet_y_lo[i][k]);
                double top    = (j < 5-1 ? result[i][j+1][k] : dirichlet_y_hi[i][k]);
                double front  = (k > 0   ? result[i][j][k-1] : dirichlet_z_lo[i][j]);
                double back   = (k < 7-1 ? result[i][j][k+1] : dirichlet_z_hi[i][j]);
                double center = result[i][j][k];

                result_laplacian[i][j][k] = (
                    ((right + left - 2 * center)) / (sidelengths[0] * sidelengths[0])
                  + ((top + bottom - 2 * center)) / (sidelengths[1] * sidelengths[1])
                  + ((front + back - 2 * center)) / (sidelengths[2] * sidelengths[2])
                );
            }
        }
    }

    // Print results

    std::cout << "Result laplacian" << std::endl;
    print_grid_3d(3, 5, 7, &result_laplacian[0][0][0]);
    std::cout << std::endl;

    std::cout
        << "Please verify that the result laplacian matches the target laplacian printed previously." << std::endl
        << "This laplacian was calculated using sidelengths "
            << sidelengths[0] << ", " << sidelengths[1] << ", and " << sidelengths[2] << " "
            << "for the x, y, and z axes respectively," << std::endl
        << "and uses the printed nonzero Dirichlet boundary values." << std::endl
        << "Please refer to the test code to verify that the laplacian is calculated as expected." << std::endl
        << std::endl
        << std::endl;
}

int main() {
    test1();
    test2();
    test3();
    test4();
}