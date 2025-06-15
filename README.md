# Linear_Algebra_Matrix

This repository provides C++17 implementations for matrix and vector calculations, covering both exact rational arithmetic and modular arithmetic over Galois Fields. No external libraries or extensions are required; all code runs with basic C++17.

---

## Features

### 1. Rational Arithmetic (`Ratio` Type)

- Exact value calculations using a custom `Ratio` type:
    - Structure: sign bit, unsigned long long denominator and numerator
    - Operators: `<, <=, >, >=, ==, !=, +, -, *, /`
    - Utility functions: `normalize()`, `Ratio_print()`, `Ratio_power()`
- **Limitations:** Square roots and non-integer exponents (e.g., `2^(1/2)`, `3^(3.4)`) are not supported. For these, use `long double`.

---

### 2. Matrix Operations (Rational Arithmetic)

- Operators: `==, !=, *, +, -`
- Functions:
    - `matrix_transpose()`
    - `RatioMat_to_LDMat()` – converts Ratio matrix to long double matrix
    - `matrix_power()`
    - `matrix_print()`
    - `matrix_inverse()`
    - `matrix_determinant()`
    - `QR_decomposition()` – Q columns are orthogonal (not normalized if using Ratio)
    - `LU_decomposition()`
    - `Eigen_Approx()` – only for symmetric, non-singular matrices (`long double`)
    - `Null_Space()`
    - `Ax=b()` – finds particular solution if multiple exist
    - `matrix_full_row_rank()` – reduces m×n matrix to r×n (or r×r if repeated with transpose)
    - `change_of_basis_P()` – computes change-of-basis matrix between two bases
- Supports vector multiplication and matrix-vector multiplication.

---

### 3. Linear Algebra over Galois Field (`Linear_Algebra_finite`)

- Modular value is user-configurable and must be a prime.
- Type: `long long` (modular arithmetic)
- Functions:
    - Modular routines: `inverse()`, `power()`, `gcd()`, `decompose()`, `divisor()`, `Extended_Euclid()`
    - Matrix operators: `==, !=, *, +, -, |` (diagonal expansion)
    - `matrix_transpose()`
    - `matrix_power()`
    - `matrix_rank()`
    - `matrix_inverse()`
    - `matrix_determinant()`
    - `Null_Space()`
    - `Ax=b()` – finds particular solution if more than one exists
    - `is_in()` – checks if vector is in column space of matrix
    - Eigenvalue decomposition:
        - `matrix_diagonalize_BF()` – brute force all values 1 to p-1
        - `matrix_diagonalize_fast()` – improved (see docs)
        - `matrix_diagonalize_henry()` – fastest (see docs)

---

### 4. Parallel Galois Field Computing (`Linear_Algebra_finite_parallel`)

- Modified version of `Linear_Algebra_finite` for parallel computing.
- Uses `omp.h` and threading to accelerate computations for large matrices (N > 200).
- Advanced workload distribution: e.g., in `matrix_diagonalize_henry`.
- **Warning:** Parallelism may cause race conditions or scheduler-related errors. For 100% accuracy, use the serial (`Linear_Algebra_finite.cpp`) version.

---

## Usage

- All functionality is accessed via the `main()` function in each file.
- Compilation does not require any non-standard headers (except `omp.h` for parallel).

Example compilation:

```
g++ -std=c++17 Linear_Algebra_Matrix.cpp -o linear_algebra_matrix
g++ -std=c++17 Linear_Algebra_finite.cpp -o linear_algebra_finite
g++ -std=c++17 -fopenmp Linear_Algebra_finite_parallel.cpp -o linear_algebra_finite_parallel

```

---

## Best Practices

- Use the `Ratio` type for exact arithmetic when integer/fraction results are required.
- For operations requiring square roots, non-integer powers, or normalization, use `long double`.
- Use the Galois field module for all modular arithmetic, ensuring your modulus is prime.
- Employ the parallel version for large matrices, but prefer the serial implementation for smaller sizes or if strict correctness is essential.

---

## Summary Table

| Module | Numeric Type | Key Features |
| --- | --- | --- |
| Linear_Algebra_Matrix | Ratio / long double | Exact matrix, vector ops; QR/LU/Eigen |
| Linear_Algebra_finite | long long (mod prime) | Modular matrix ops; null space; eigenvalues |
| Linear_Algebra_finite_parallel | long long (mod prime) | Parallel matrix ops for large N |

---

## Documentation

- Every function is documented with comments.
- For deeper mathematical and algorithmic explanations, see the inline documentation and (where referenced) supplementary documents in the repository.

---

## Acknowledgement

Good luck with your computations and mathematical research!
