#include "fd_gauge_correction.h"
#include <iostream>
#include <string>

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
    print_title("Test 2: Poisson solve");

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


    double sidelengths[3] = { 1, 2, 4 };

    std::cout << "Target laplacian" << std::endl;
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
        << "This laplacian was calculated uses sidelengths "
            << sidelengths[0] << ", " << sidelengths[1] << ", and " << sidelengths[2] << " "
            << "for the x, y, and z axes respectively." << std::endl
        << "Please refer to the test code to verify that the laplacian is calculated as expected." << std::endl;
}

int main() {
    test1();
    test2();
}