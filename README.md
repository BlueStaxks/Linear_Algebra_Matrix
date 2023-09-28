# Linear_Algebra_Matrix
 


This code is for Matrix and Vector calculation.

This code has "Ratio" type. It has sign bit, unsigned long long denominator and numerator.
+,-,*,/ operation is defined.

You can use this Ratio type to calculate exact value (bounded by unsigned long long).

Unfortunately, Ratio type cannot handle square root. For example, 2^(1/2) or 3^(3.4).
This is very unfortunate since QR decomposition and Eigenvalues usually require square root.
For these applications, you need to use long double type.



This code includes following functions.

For Ratio type :

normalize()
<
<=
>
>=
==
!=
Ratio_print()
*
/
+
-
Ratio_power()

-----------------------

For Matrix : 

==
!=
*
+
-
matrix_transpose()
RatioMat_to_LDMat()      // This function converts Ratio type Matrix to long double type Matrix
matrix_power()
matrix_print()
matrix_inverse()
matrix_determinant()
QR_decomposition()       // If you use Ratio type Matrix, Q's col vectors are only orthogonal. Not unit vectors. But QR is still A.
LU_decomposition()
Eigen_Approx()           // only Symmetric Matrix can be used.



Vector's multiplication and Matrix_vector multiplication are also defined.
