#ifndef PROJECT1_DIRECT_H
#define PROJECT1_DIRECT_H

// Include header file
#include "dense_linear_algebra.h"

// Import necessary libraries
#include <stdexcept> // Library to manage exceptions
#include <cmath>     // Library for maths operations


// Use previously declared namespace for linear algebra to add new ojects there
namespace BasicDenseLinearAlgebra {

// Base class for a linear colver
class LinearSolver {
public:
    // Solves the linear system Ax = b and returns vector x as a solution
    // Function is pure virtual
    virtual DoubleVector solve(const SquareDoubleMatrix& A, const DoubleVector& b) const = 0;
    
    // Destructor to avoid memory leakages
    virtual ~LinearSolver() {}
};


// Class implementing LU decomposition, taken from dense_linear_algebra.h header file
class LUSolver : public LinearSolver {
public:
    // Function that solves the linear system using LU decomposition
    // Overrides base class function
    DoubleVector solve(const SquareDoubleMatrix& A, const DoubleVector& b) const override {
        // Create LU Solver object
        LULinearSolver luSolver;

        // Solve the system and return the solution vector
        return luSolver.lu_solve(A, b);
    }
};

// Class implementing Cholesky decomposition
class CholeskySolver : public LinearSolver {
public:
    // Overrides base class function
    DoubleVector solve(const SquareDoubleMatrix& A, const DoubleVector& b) const override {
        // Check if the matrix is symmetric using method from the private
        // We don't need to check if the matrix is square as it is intended to be square
        // Because SquareDoubleMatrix class enforces this requirement
        if (!isSymmetric(A)) {
            throw LinearSolverError("Matrix is not symmetric.");
        }

        // Get matrix size
        unsigned n = A.n();

        // Initialize a vector to store lower triangular matrix L
        std::vector<double> L(n * n, 0.0);


        // Loop to perform Cholesky decomposition
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j <= i; ++j) {
                // Start with the first element of the matrix
                double sum = A(i, j);

                // Calculate the sum
                for (unsigned k = 0; k < j; ++k) {
                    sum -= L[i * n + k] * L[j * n + k];
                }
                // Check if the matrix is positive definite
                if (i == j) {
                    if (sum <= 0.0) {
                        throw LinearSolverError("Matrix is not positive definite.");
                    }
                    // Calculate the square root of the diagonal
                    L[i * n + j] = std::sqrt(sum);
                
                // Calculate off-diagonal element
                } else {
                    L[i * n + j] = sum / L[j * n + j];
                }
            }
        }
        
        // Forward substitution to solve L * y = b
        // Vector to store solution values 
        DoubleVector y(n);
        for (unsigned i = 0; i < n; ++i) {
            // Start with right-hand side vector
            double sum = b[i];

            // Substract known values of vector y
            for (unsigned j = 0; j < i; ++j) {
                sum -= L[i * n + j] * y[j];
            }

            // Divide by the diagonal element of L 
            y[i] = sum / L[i * n + i];
        }

        // Backward substitution to solve L^T * x = y 
        // Initialiaze solution vector
        DoubleVector x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[i];

            // Substract known values of vector x
            for (unsigned j = i + 1; j < n; ++j) {
                sum -= L[j * n + i] * x[j];
            }
            // Divide by the diagonal element of L 
            x[i] = sum / L[i * n + i];
        }

        // Return solution vector
        return x;
    }

private:
    // Function to check if the matrix is symmetric
    // If all the values A(i, j) == A(j, i), matrix is symmetric
    bool isSymmetric(const SquareDoubleMatrix& A) const {
        
        // Get matrix size
        unsigned n = A.n();

        // Loop to check condition stated above
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < i; ++j) {
                if (A(i, j) != A(j, i)) {
                    // Return False if condition is not satisfied
                    return false;
                }
            }
        }
        // Return True if condition is satisfied
        return true;
    }

};

} // End of the namespace

#endif // PROJECT1_DIRECT_H