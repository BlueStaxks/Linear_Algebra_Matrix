# Linear_Algebra_Matrix
 


This code is for Matrix and Vector calculation.

This code has "Ratio" type. It has sign bit, unsigned long long denominator and numerator.




You can use this Ratio type to calculate exact value (bounded by unsigned long long).

Unfortunately, Ratio type cannot handle square root. For example, 2^(1/2) or 3^(3.4).

This is very unfortunate since QR decomposition and Eigenvalues usually require square root.

For these applications, you need to use long double type.



----------------------
-----------------------

This code includes following functions.

For Ratio type :

<, <=, >, >=, ==, !=, *, /, +, -

normalize(), Ratio_print(), Ratio_power()

-----------------------

For Matrix : 

==, !=, *, +, -


matrix_transpose()

RatioMat_to_LDMat()      // This function converts Ratio type Matrix to long double type Matrix

matrix_power(), matrix_print()

matrix_inverse()

matrix_determinant()

QR_decomposition()       // If you use Ratio type Matrix, Q's col vectors are only orthogonal. Not unit vectors. But QR is still A.

LU_decomposition()

Eigen_Approx()           // only Symmetric Matrix can be used(Not Singular).

Null_Space()

Ax=b()                   // This function finds x. if there is more than one solution, this function will return particular solution.

matrix_full_row_rank()   // Makes m by n matrix to r by n matrix. (Do it and do it with transpose matrix, then it will become r by r invertible matrix or scalar)

change_of_basis_P()      // This function uses 2 basis and make change-of-basis (B to C) matrix P.

--------------------------

Vector's multiplication and Matrix_vector multiplication are also defined.

--------------------------
------------------------
---------------------------


# Linear_Algebra_Matrix(Galois Field)
 

Linear_Algebra_finite is about Galois field.

You can set the modular value yourself. it must be a prime number.

The code will run Initiation() with your modular number.

only long long type is available.



----------------------
-----------------------

This code includes following functions.

For long long type :

inverse(), power(), gcd(), decompose(), divisior(), Extended_Euclid()

-----------------------

For Matrix : 

==, !=, *, +, -, |(diagonal expansion)

matrix_transpose()

matrix_power(), matrix_print()

matrix_rank()

matrix_inverse()

matrix_determinant()

Null_Space()

Ax=b()                     // This function finds x. if there is more than one solution, this function will return particular solution.

matrix_diagonalize_BF()    // This function do brute force way to find all eigenvalues. Trying 1 to p-1.

matrix_diagonalize_fast()  // This function is a faster diagonalize function. check the document about Galois Field Matrix in this repository. (The document is not yet fully written)

matrix_diagonalize_henry() // This is the fasteset way I found to diagonalize. check the document for mathematical detail. 

---------------------

Vector's multiplication and Matrix_vector multiplication are also defined.

Use main() to do your calculations!!

Good Luck
