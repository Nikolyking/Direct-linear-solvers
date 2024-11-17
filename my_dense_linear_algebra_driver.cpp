#include <iostream>
#include "project1_direct.h"


// Test functions to check if the answer of LU solver is correct
void testLUSolver() {
    // Define a 3x3 matrix A
    BasicDenseLinearAlgebra::SquareDoubleMatrix A(3);
    A(0, 0) = 2.0; A(0, 1) = 2.0; A(0, 2) = 0.0;
    A(1, 0) = 2.0; A(1, 1) = 2.0; A(1, 2) = 0.0;
    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 2.0;

    // Define the right-hand side vector b
    BasicDenseLinearAlgebra::DoubleVector b(3);
    b[0] = 1.0;
    b[1] = 1.0;
    b[2] = 1.0;

    // Expected solution vector x
    BasicDenseLinearAlgebra::DoubleVector expected_x(3);
    expected_x[0] = 1.0;
    expected_x[1] = 1.0;
    expected_x[2] = 1.0;

    // Solve Ax = b using LUSolver
    BasicDenseLinearAlgebra::LUSolver luSolver;
    BasicDenseLinearAlgebra::DoubleVector x = luSolver.solve(A, b);

    // Output and check the results
    std::cout << "LUSolver solution:" << std::endl;
    for (unsigned i = 0; i < x.n(); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
        assert(std::abs(x[i] - expected_x[i]) < 1e-6); // Check with a small tolerance
    }
}


// Test functions to check if the answer of Cholesky solver is correct
void testCholeskySolver() {
    // Define a 3x3 symmetric positive-definite matrix A
    BasicDenseLinearAlgebra::SquareDoubleMatrix A(3);
    A(0, 0) = 2.0; A(0, 1) = 2.0; A(0, 2) = 0.0;
    A(1, 0) = 2.0; A(1, 1) = 2.0; A(1, 2) = 0.0;
    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 2.0;

    // Define the right-hand side vector b
    BasicDenseLinearAlgebra::DoubleVector b(3);
    b[0] = 1.0;
    b[1] = 1.0;
    b[2] = 1.0;

    // Expected solution vector x
    BasicDenseLinearAlgebra::DoubleVector expected_x(3);
    expected_x[0] = 1.0;
    expected_x[1] = 1.0;
    expected_x[2] = 1.0;

    // Solve Ax = b using CholeskySolver
    BasicDenseLinearAlgebra::CholeskySolver choleskySolver;
    BasicDenseLinearAlgebra::DoubleVector x = choleskySolver.solve(A, b);

    // Output the result
    std::cout << "Chelesky solution:" << std::endl;
    for (unsigned i = 0; i < x.n(); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
}

int main() {
    // testLUSolver();
    testCholeskySolver();
    
    return 0;
}