#ifndef DENSE_LINEAR_ALGEBRA_H
#define DENSE_LINEAR_ALGEBRA_H


// C++ includes
#include<vector>
#include<string>
#include<fstream>
#include<cmath>
#include<cassert>
#include<sstream>
#include<chrono>


// ############################################################
/// Helper file with functions/classes for basic linear algebra
// ############################################################



//============================================================
/// Namespace for basic dense linear algebra
//============================================================
namespace BasicDenseLinearAlgebra
{
 
//===================================================
/// Class for a dense general (not necessarily
/// square) matrix of doubles
//===================================================
 class DoubleMatrix
 {
 
 public:

  /// Constructor: Pass size (n rows, m columns)
  DoubleMatrix(const unsigned& n, const unsigned& m) : N(n), M(m)
   {
    // Flat packed storage: a(i,j) = a_flat_packed(i*M+j) (row by row)
    Matrix_storage.resize(N*M,0.0);
   }
  
  /// Number of rows
  unsigned n() const
   {
    return N;
   }
  
  /// Number of columns
  unsigned m() const
   {
    return M;
   }
 
  /// Const access to (i,j)-th entry
  double operator()(const unsigned& i, const unsigned& j) const
   {
#ifdef RANGE_CHECKING
    assert(i<N);
    assert(j<M);
#endif
    return Matrix_storage[i*M+j];
   }

  /// Read/write access to (i,j)-th entry
  double& operator()(const unsigned& i, const unsigned& j) 
   {
#ifdef RANGE_CHECKING
    assert(i<N);
    assert(j<M);
#endif
    return Matrix_storage[i*M+j];
   }

  /// Output to std::cout
  void output() const
   {
    output(std::cout);
   }

  /// Output to file (specify filename)
  void output(std::string filename) const
   {
    std::ofstream outfile(filename.c_str());
    output(outfile);
    outfile.close();
   }

  
  /// Output to file (specify stream)
  void output(std::ostream& outfile) const
   {
    for (unsigned i=0;i<N;i++)
     {
      for (unsigned j=0;j<M;j++)
       {
        outfile << i << " "
                << j << " "
                << Matrix_storage[i*M+j] << " "
                << std::endl;
       }
     }
   }
  
  /// Input from file (specify filename)
  void read(std::string filename)
   {
    std::ifstream infile(filename.c_str());
    read(infile);
    infile.close();
   }

  
  /// Input from file (specify stream)
  void read(std::ifstream& infile)
   {
    unsigned i_read=0;
    unsigned j_read=0;
    for (unsigned i=0;i<N;i++)
     {
      for (unsigned j=0;j<M;j++)
       {
        infile >> i_read;        
        if (i!=i_read)
         {       
          std::stringstream str;
          str << "\n\nERROR: Row index in matrix is i = "
              << i <<  " but data is provided for " << i_read
              << std::endl;
          throw std::runtime_error(str.str().c_str());
         }
        infile >> j_read;        
        if (j!=j_read)
         {       
          std::stringstream str;
          str << "\n\nERROR: Column index in matrix is j = "
              << i <<  " but data is provided for " << j_read
              << std::endl;
          throw std::runtime_error(str.str().c_str());
         } 
        infile >> Matrix_storage[i*M+j];
       }
     }
   }

    
 private:
 
  /// Number of rows
  unsigned N;
  
  /// Number of columns
  unsigned M;

  /// Entries are flat packed, row by row:
  /// a(i,j) = a_flat_packed(i*M+j) (row by row)
  std::vector<double> Matrix_storage;

  
 };


 
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

 
//===================================================
/// Class for a dense square matrix of doubles
//===================================================
 class SquareDoubleMatrix : public DoubleMatrix
 {

 public:

  /// Constructor: Pass size
  SquareDoubleMatrix(const unsigned& n) : DoubleMatrix(n,n)
   {}
  
 };
  


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

 
//===================================================
/// Class for a vector of doubles
//===================================================
 class DoubleVector
 {
 
 public:
 
  /// Constructor: Pass size
  DoubleVector(const unsigned& n=0) : N(n)
   {
    // Resize storage and initialise entries to zero
    Vector_storage.resize(N,0.0);
   }

  
  /// Size of vector
  unsigned n() const
   {
    return N;
   }
 
  /// Const access to i-th entry
  double operator[](const unsigned& i) const
   {
#ifdef RANGE_CHECKING
    assert(i<N);
#endif
    return Vector_storage[i];
   }

  /// Read/write access to i-th entry
  double& operator[](const unsigned& i)
   {
#ifdef RANGE_CHECKING
    assert(i<N);
#endif
    return Vector_storage[i];
   }

 /// Resize (and zero the entries)
 void resize(const unsigned& n)
  {
   N=n;
   Vector_storage.resize(n,0.0);
  }


  /// Output to std::cout
  void output() const
   {
    output(std::cout);
   }
  
  /// Output (specify filename)
  void output(std::string filename) const
   {
    std::ofstream outfile(filename.c_str());
    output(outfile);
    outfile.close();
   }
  
  /// Output (specify stream)
  void output(std::ostream& outfile) const
   {
    for (unsigned i=0;i<N;i++)
     {
      outfile << i << " "
              << Vector_storage[i] << " " 
              << std::endl;
     }
   }
 
 
  /// Input from file (specify filename)
  void read(std::string filename)
   {
    std::ifstream infile(filename.c_str());
    read(infile);
    infile.close();
   }

  
  /// Input from file (specify stream)
  void read(std::ifstream& infile)
   {
    unsigned i_read=0;
    for (unsigned i=0;i<N;i++)
     {
      infile >> i_read;        
      if (i!=i_read)
       {       
        std::stringstream str;
        str << "\n\nERROR: Row index in matrix is i = "
            << i <<  " but data is provided for " << i_read
            << std::endl;
        throw std::runtime_error(str.str().c_str());
       }
      infile >> Vector_storage[i];
     }
   }
  
 protected:

  /// Size of the vector
  // Note: yes, this duplicates information
  // but is quicker than calling Vector_storage.size().
  // It's safe to do this since the user can't over-write
  // N from the outside so the two stay in sync.
  unsigned N;
 
  /// Entries are stored internally as a std::vector of doubles
  std::vector<double> Vector_storage;
 
 };




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


 
//===================================================
/// Error for linear solvers
//===================================================
 class LinearSolverError : public std::runtime_error
 {
  
 public:

  /// Issue runtime error, outputting generic message
  LinearSolverError() :
   runtime_error("Error in linear solver!")
   {}

  /// Runtime error with context specific message
  LinearSolverError(std::string msg) :
   runtime_error(msg.c_str())
   {}


 };


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

 
 
//=============================================================================
/// Dense LU decomposition-based solver
//============================================================================
 class LULinearSolver
 {

 public:

  /// Constructor
  LULinearSolver() {}

  /// Destructor
  ~LULinearSolver() {}


  /// Do the linear solve: Takes matrix and rhs
  /// vector and returns the solution of the linear system.
  /// (Not const because it updates some internal storage for
  /// LU factors etc.)
  DoubleVector lu_solve(const SquareDoubleMatrix&  matrix,
                        const DoubleVector& rhs); 

 private:

  /// Perform the LU decomposition of the matrix
  void factorise(const SquareDoubleMatrix& matrix);
 
  /// Do the backsubstitution step to solve the system LU result = rhs
  DoubleVector backsub(const DoubleVector& rhs);

  /// Storage for the index of permutations in the LU solve
  /// (used to handle pivoting)
  std::vector<unsigned> Index;
 
  /// Storage for the LU decomposition (flat-packed into nxn vector)
  std::vector<double> LU_factors;
 
 };




//=============================================================================
/// Linear solver: Takes matrix and rhs
/// vector and returns the solution of the linear system.
//============================================================================
 DoubleVector LULinearSolver::lu_solve(const SquareDoubleMatrix& matrix,
                                       const DoubleVector& rhs)
 {
  // factorise
  factorise(matrix);
  
  // Get result via backsubstitution
  DoubleVector result=backsub(rhs);

  return result;
 }




//=============================================================================
/// LU decompose the matrix.
//=============================================================================
 void LULinearSolver::factorise(const SquareDoubleMatrix& matrix)
 {
  // Set the number of unknowns
  const unsigned n = matrix.n();
   
  // Allocate storage for the LU factors and the permutation index
  // set entries to zero.
  LU_factors.resize(n*n,0.0);
  Index.resize(n,0);
 
  // Now we know that memory has been allocated, copy over
  // the matrix values
  unsigned count = 0;
  for (unsigned i = 0; i < n; i++)
   {
    for (unsigned j = 0; j < n; j++)
     {
      LU_factors[count] = matrix(i, j);
      ++count;
     }
   }
 
  // Loop over columns
  for (unsigned j = 0; j < n; j++)
   {
    // Initialise imax, the row with the largest entry in the present column
    unsigned imax = 0;

    // Do rows up to diagonal
    for (unsigned i = 0; i < j; i++)
     {
      double sum = LU_factors[n * i + j];
      for (unsigned k = 0; k < i; k++)
       {
        sum -= LU_factors[n * i + k] * LU_factors[n * k + j];
       }
      LU_factors[n * i + j] = sum;
     }
   
    // Initialise search for largest pivot element
    double largest_entry = 0.0;

    // Do rows below diagonal -- here we still have to pivot!
    for (unsigned i = j; i < n; i++)
     {
      double sum = LU_factors[n * i + j];
      for (unsigned k = 0; k < j; k++)
       {
        sum -= LU_factors[n * i + k] * LU_factors[n * k + j];
       }
      LU_factors[n * i + j] = sum;
      
      // New largest entry found in a row below the diagonal?
      double tmp = std::fabs(sum);
      if (tmp >= largest_entry)
       {
        largest_entry = tmp;
        imax = i;
       }
     }
   
    // Test to see if we need to interchange rows; if so, do it!
    if (j != imax)
     {
      for (unsigned k = 0; k < n; k++)
       {
        double tmp = LU_factors[n * imax + k];
        LU_factors[n * imax + k] = LU_factors[n * j + k];
        LU_factors[n * j + k] = tmp;
       }
          
     }
   
    // Record the index (renumbering rows of the orignal linear
    // system to reflect pivoting)
    Index[j] = imax;
   
    // Divide by pivot element
    if (j != n - 1)
     {
      double pivot= LU_factors[n * j + j];
      if (pivot==0.0)
       {
        std::string error_message=
         "Singular matrix: zero pivot in row "+std::to_string(j);
        throw LinearSolverError(error_message.c_str());
       }
      double tmp = 1.0 / pivot;
      for (unsigned i = j + 1; i < n; i++)
       {
        LU_factors[n * i + j] *= tmp;
       }
     }
   
   } // End of loop over columns
 
 
 }
 

//=============================================================================
/// Do the backsubstitution for the DenseLU solver.
//=============================================================================
 DoubleVector LULinearSolver::backsub(const DoubleVector& rhs)
 {
 
  // Initially copy the rhs vector into the result vector
  const unsigned n = rhs.n();
  DoubleVector result(n);
  for (unsigned i = 0; i < n; ++i)
   {
    result[i] = rhs[i];
   }
 
  // Loop over all rows for forward substitution
  unsigned k = 0;
  for (unsigned i = 0; i < n; i++)
   {
    unsigned ip = Index[i];
    double sum = result[ip];
    result[ip] = result[i];
    if (k != 0)
     {
      for (unsigned j = k - 1; j < i; j++)
       {
        sum -= LU_factors[n * i + j] * result[j];
       }
     }
    else if (sum != 0.0)
     {
      k = i + 1;
     }
    result[i] = sum;
   }
 
  // Now do the back substitution
  // Note: this has to be an int to avoid wrapping around!
  for (int i = n - 1; i >= 0; i--)
   {
    double sum = result[i];
    for (unsigned j = i + 1; j < n; j++)
     {
      sum -= LU_factors[n * i + j] * result[j];
     }
    result[i] = sum / LU_factors[n * i + i];
   }

  return result;
 }




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



//==========================================================================
/// Helper function to get the max. error of the solution of Ax=b, defined
/// as max_i |(A_ij x_j - b_i)|
//========================================================================== 
 double max_error(const SquareDoubleMatrix& matrix,
                  const DoubleVector& rhs,
                  const DoubleVector& soln)
 {
  unsigned n=rhs.n();
  double max_error=0.0;
  for (unsigned i=0;i<n;i++)
   {
    double error=rhs[i];
    for (unsigned j=0;j<n;j++)
     {
      error-=matrix(i,j)*soln[j];
     }
    if (fabs(error)>max_error) max_error=error;
   }
  return max_error;
 }
 
 
} // end of namespace
  
  

#endif