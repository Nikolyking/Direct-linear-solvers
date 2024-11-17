#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include "project1_direct.h"

using namespace BasicDenseLinearAlgebra;

/////////////////////////////////////////////////////
/////////// Functions for multiple tests ////////////
/////////////////////////////////////////////////////

// Function to generate a symmetric positive definite matrix
SquareDoubleMatrix generateSPDMatrix(unsigned n) {
    SquareDoubleMatrix A(n);

    // Loop over rows and columns
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j <= i; ++j) {
            
            // Fill matrix with normalized random numbers between 0 and 1
            double value = rand() / (double)RAND_MAX;
            A(i, j) = value;
            A(j, i) = value;
        }

        // Add some number to the diagonal to ensure dominance
        A(i, i) += n;
    }
    return A;
}

// Function to generate a general square matrix
SquareDoubleMatrix generateGeneralMatrix(unsigned n) {
    SquareDoubleMatrix A(n);

    // Loop over rows and columns
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {

            // Fill matrix with normalized random numbers between 0 and 1
            A(i, j) = rand() / (double)RAND_MAX;
        }
    }
    return A;
}

int main() {
    // Initialize vector with matrix sizes
    std::vector<unsigned> sizes = {100, 200, 300, 400, 500, 1000, 2000, 3000};

    // Initialize vectors for consumed time computations 
    std::vector<double> lu_times, chol_times;

    // Loop over vector
    for (const unsigned& n: sizes) {
        // Generate matrixes
        SquareDoubleMatrix A_lu = generateGeneralMatrix(n);
        SquareDoubleMatrix A_chol = generateSPDMatrix(n);

        // Generate the right hand-side vector with values from 0 to 1
        DoubleVector b(n);
        for (unsigned i = 0; i < n; ++i) {
            b[i] = rand() / (double)RAND_MAX;
        }

        // Creare LUSolver object
        LUSolver lu_solver;

        // Initialize chrono and perform computation
        auto start_lu = std::chrono::high_resolution_clock::now();
        DoubleVector x_lu = lu_solver.solve(A_lu, b);
        auto end_lu = std::chrono::high_resolution_clock::now();

        // Write calculation time
        std::chrono::duration<double> elapsed_lu = end_lu - start_lu;
        lu_times.push_back(elapsed_lu.count());

        // Creare CholeskySolver object
        CholeskySolver chol_solver;

        // Initialize chrono and perform computation
        auto start_chol = std::chrono::high_resolution_clock::now();
        DoubleVector x_chol = chol_solver.solve(A_chol, b);
        auto end_chol = std::chrono::high_resolution_clock::now();

        // Write calculation time
        std::chrono::duration<double> elapsed_chol = end_chol - start_chol;
        chol_times.push_back(elapsed_chol.count());
    }

    // Output results to a file
    std::ofstream outfile("solve_times.txt");
    outfile << "n lu_time chol_time\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        outfile << sizes[i] << " " << lu_times[i] << " " << chol_times[i] << "\n";
    }
    outfile.close();

    return 0;
}