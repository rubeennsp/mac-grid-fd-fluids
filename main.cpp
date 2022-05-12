#include "fd_gauge_correction.h"
#include <iostream>

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

    double grid[3][5][7] = {
        {
            {111, 112, 113, 114, 115, 116, 117, }, 
            {121, 122, 123, 124, 125, 126, 127, }, 
            {131, 132, 133, 134, 135, 136, 137, }, 
            {141, 142, 143, 144, 145, 146, 147, }, 
            {151, 152, 153, 154, 155, 156, 157, },
        },
        {
            {211, 212, 213, 214, 215, 216, 217, }, 
            {221, 222, 223, 224, 225, 226, 227, }, 
            {231, 232, 233, 234, 235, 236, 237, }, 
            {241, 242, 243, 244, 245, 246, 247, }, 
            {251, 252, 253, 254, 255, 256, 257, },
        },
        {
            {311, 312, 313, 314, 315, 316, 317, }, 
            {321, 322, 323, 324, 325, 326, 327, }, 
            {331, 332, 333, 334, 335, 336, 337, }, 
            {341, 342, 343, 344, 345, 346, 347, }, 
            {351, 352, 353, 354, 355, 356, 357, },
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

int main() {
    test1();
}