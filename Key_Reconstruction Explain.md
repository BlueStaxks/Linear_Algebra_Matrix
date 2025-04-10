### Introduce

The code is designed to detect fake shares and successfully reconstruct the secret key K.

This detection algorithm is highly robust. It can reconstruct the key even with 11 valid shares and 1000 fake shares (t = 10, s = 1001, prime = 524287).

This method requires an (s+1) × m matrix, where s is the maximum number of fake shares, and m is the total number of shares (valid and fake).

The receiver uses the matrix to identify valid shares and reconstructs K using only these identified valid shares.

If the number of real shares is fewer than t+1, the system cannot reconstruct K because the matrix does not contain information about the y-coordinates, thus ensuring security.

*Note: This method is not guaranteed to always succeed, but it has an extremely high probability of success.

---

### Method:

1. Side A initiates by generating n Shamir shares, including both valid and fake shares, in a Galois Field (GF) with prime p. Any subset of t+1 valid shares can reconstruct the secret key.

2. Side A creates an s × t matrix in GF. Each row corresponds to a polynomial of degree at most t-1, with s representing the maximum number of fake shares.

3. Side A calculates a vector for each share as follows:
   - 3.1. If a share is valid, the vector is defined as r = ⟨x, p₁(x), p₂(x), …, pₛ(x)⟩, where x is the x-coordinate of the share.
   - 3.2. If a share is fake, the vector is defined as r = ⟨x, r₁, r₂, …, rₛ⟩, where x is the x-coordinate, and each rᵢ is a random number.

4. Side A organizes these vectors as columns of matrix A.

5. Side A transmits matrix A along with the n Shamir shares to Side B.

---

Because Side A selects the shares' x-coordinates randomly, columns of matrix A cannot be distinguished as valid or fake when fewer than t+1 valid shares are present. If n > t, reconstruction attempts may fail due to insufficient valid shares. Side A deliberately introduces ambiguity about the real share count by injecting fake shares, making it impossible for Side B to ascertain the exact number of real shares unless more than t valid shares exist.

---

6. Side B receives matrix A and the n Shamir shares. Both sides know values s and t.

7. Side B expands matrix A into matrix M:
   - 7.1. Extract the first row of A and name it vector v.
   - 7.2. Append t rows on top of A, where the ith row is v^i, with i ranging from 0 to t-1. (v^i means raising each element of v to the power i.)
   - 7.2. Name it M.

8. Side B computes the Null Space of matrix M. A trivial Null Space indicates insufficient valid shares.

9. Side B identifies indices of columns corresponding to non-zero entries in vectors of the Null Space.

10. Side B extracts the first elements of the identified columns from matrix A; these represent the x-coordinates of valid shares.

11. Side B reconstructs the secret key K using these valid shares.

---

### Analysis of Gaussian Elimination and Probability of Failure

Since matrix M incorporates the first t rows generated from vector v, Gaussian elimination can completely remove the last s rows derived from genuine polynomial shares. This occurs because linear combinations of these first t rows can represent any polynomial of degree at most t-1. In contrast, columns derived from fake shares cannot be fully eliminated.

If we assume an ideal scenario with all fake columns placed on the left, the resulting matrix after clean elimination would look like:

C = [X; R 0],

where R represents the bottom-left random part submatrix.

Of course, such a perfect elimination is not possible for side B because Side B does not know the key matrix created by Side A. Nevertheless, since Gaussian elimination does not change the NullSpace of the matrix, the NullSpace computed by Side B is identical to the null space of matrix C.

It can be easily seen that the NullSpace of C has zero entries corresponding to fake columns if these columns are random vectors. This clarifies why at least t+1 real columns are required for the existence of a non-trivial NullSpace.

If only t real columns exist, the problem reduces to examining the top-right t by t submatrix X from C. Matrix X is invertible because its rows are constructed to be linearly independent.

Consequently, after Side B computes the NullSpace, a trivial NullSpace implies fewer than t+1 valid shares, as at least t+1 valid shares are necessary for a non-trivial NullSpace.



The method may fail if a randomly generated fake share coincidentally matches a valid polynomial column. The probability of this occurring for a single fake column is (1/p)^s, given that s random numbers must match the actual polynomial values p₁(x), ..., pₛ(x).

Thus, the overall failure probability is calculated as:

Failure Probability = 1 - (Probability all fake columns differ from valid columns)  
= 1 - (1 - (1/p)^s)^(number of fake shares)

Even with s fake shares, the probability of failure rapidly approaches zero:
- For example, when p = 524287 and s = 1, the failure probability ≈ 0.000002.
- For s ≥ 2, the failure probability is virtually negligible.

Even though the probability of failure is negligible, this issue can be avoided by verifying the randomly chosen vector and replacing it if it coincides with an actual vector from s polynomials.
