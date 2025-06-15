Linear Algebra Library
This repository offers a comprehensive suite of matrix and vector operations in C++17, supporting both exact rational arithmetic and computations over finite (Galois) fields. The design emphasizes clarity, correctness, and academic rigor, suitable for teaching, research, and practical engineering.

1. Exact Arithmetic (Ratio Type)
Implements exact rational numbers using:

struct Ratio {
    bool sign;
    unsigned long long num;
    unsigned long long den;
};
Features:
Arithmetic & comparisons: <, <=, >, >=, ==, !=, +, -, *, /

Utility routines:

normalize() — reduce to lowest terms

Ratio_print() — display ratio as readable output

Ratio_power() — integer exponentiation

⚠️ Limitations: Does not support irrational or non-integer powers (e.g., √2 or 3⁰·⁴); convert to long double when needed.

2. Matrix Operations with Ratio
Supports exact integer-based matrices and vectors:

Operators: ==, !=, +, -, * (matrix multiplication)

Standard functions:

matrix_transpose()

matrix_power()

matrix_determinant()

matrix_inverse()

matrix_print()

Decompositions & solvers:

LU_decomposition()

QR_decomposition()

Q has orthogonal columns (not normalized due to exact arithmetic), but A = Q·R holds.

Eigen_Approx() — approximates eigenvalues (only for symmetric, non-singular matrices), must convert to long double

Null_Space() — basis for null space

Ax=b() — solve linear system; returns one particular solution if under‑determined

matrix_full_row_rank() — reduce an m×n matrix to full row rank

change_of_basis_P() — compute basis transition matrix between two given bases

Additionally, vector–vector and matrix–vector operations are fully implemented.

3. Galois Field Arithmetic
Supports arithmetic over ℤₚ (prime modulus):

Initialization: Initiation(p) where p is a prime

Modular arithmetic utilities:
inverse(), power(), gcd(), Extended_Euclid(), decompose(), divisor()

Matrix operations over ℤₚ:
Operators: ==, !=, +, -, * (with optional diagonal expansion via |)

Standard functions:

matrix_transpose()

matrix_power()

matrix_rank()

matrix_determinant()

matrix_inverse()

Null_Space(), Ax=b() — solves, returns a particular solution if not unique

is_in() — checks if a vector belongs to the column space of a given matrix

Eigen-decomposition methods:
matrix_diagonalize_BF() — brute-force search through {1, … p−1}

matrix_diagonalize_fast() — more efficient algorithm (documented separately)

matrix_diagonalize_henry() — fastest method available (see documentation)

4. Parallel Finite-Field Variant
A thread‑parallel version using OpenMP and std::thread, optimized for matrices with N ≥ 200. It accelerates core operations (e.g., matrix_diagonalize_henry()), but may introduce non-determinism or race conditions on smaller sizes. For guaranteed correctness, prefer the serial finite‑field implementation.

5. Conversion & Mixed-Precision Usage
Direct conversion: RatioMat_to_LDMat() — from rational matrix to long double matrix

Use exact arithmetic for tasks requiring exactness

Use long double for operations involving irrational numbers (e.g., normalization in QR, eigenvalue approximation)

6. Compilation & Usage
All code is compatible with standard C++17, without external libraries:


g++ -std=c++17 Linear_Algebra_Matrix.cpp -o linear_exact
g++ -std=c++17 -fopenmp Linear_Algebra_finite_parallel.cpp -o linear_finite_parallel
linear_exact – rational-based exact arithmetic

linear_finite_parallel – parallel finite-field computations (requires -fopenmp)

Call main() in each file to use or test routines directly.

7. Best Practices
For exactness (no rounding), use Ratio and rational-matrix routines.

Where irrational numbers or approximate eigenvalues are needed, convert via RatioMat_to_LDMat() before using floating-point methods.

Use finite-field routines for modular arithmetic in ℤₚ.

For parallel matrix methods, ensure N ≥ 200 and be aware of possible concurrency variability.

8. Summary of Contents
Rational (Ratio) module:
Ratio type + arithmetic/comparison operators

normalize(), Ratio_print(), Ratio_power()

Rational Matrix module:
Matrix operators and matrix_transpose(), matrix_power(), matrix_determinant(), matrix_inverse(), matrix_print()

LU_decomposition(), QR_decomposition(), Eigen_Approx(), Null_Space(), Ax=b(), matrix_full_row_rank(), change_of_basis_P()

Finite-Field (long long mod p) module:
Initiation(p), modular arithmetic (inverse(), power(), gcd(), Extended_Euclid(), decompose(), divisor())

Matrix arithmetic, transpose, power, rank, determinant, inverse, null space, Ax=b(), is_in()

Eigen decomposition: matrix_diagonalize_BF(), matrix_diagonalize_fast(), matrix_diagonalize_henry()

Parallel Finite-Field (Linear_Algebra_finite_parallel):
Parallel acceleration via OpenMP or threads

Includes threaded eigen routines; best for large matrices

Use serial version for small sizes or precise results.

9. Glossary
Ratio: Exact rational arithmetic

Exact Matrix: Dense matrices over Ratio

Finite-Field Matrix: Dense matrices over Zₚ

Parallel Finite-Field: Threaded version for large matrices

10. Reference & Documentation
Please refer to code comments and the repository for details on algorithms, workflows, and examples.

Good luck with your computations, proofs, and implementations!
