### INTRODUCTION

The code is developed to detect and filter out fake shares and then reconstruct the secret key K using only the valid (real) shares. The algorithm is robust—it has been demonstrated in simulation that it can correctly reconstruct the key even when there are far more fake shares than valid ones (for example, with t = 10, s = 1001, using a prime of 524287, it can reconstruct the key with 11 valid shares amid 1000 fake ones).

In addition, the protocol is designed so that if fewer than t+1 valid shares are present, not only is reconstruction of K impossible but the receiver also cannot determine the exact number of real shares. This intentional ambiguity about the count of valid shares in low-sampling scenarios contributes to the overall security of the system by preventing an adversary (or even Side B, the receiver, prior to sufficient evidence) from ascertaining how many genuine shares there are.

The method requires that the receiver obtains a matrix of dimensions (s+1) × m, where s is the maximum allowed number of fake shares and m is the overall number of shares (including both valid and fake). Using this matrix, the receiver identifies which shares are genuine and then reconstructs K using just those shares. An important security feature is that if there are fewer than t+1 valid shares, the matrix lacks the necessary information in the y-coordinates for reconstruction, preventing key recovery.

A note is made that while the probability of failure is extremely low.

----

### METHOD

Step by step, the protocol operates as follows:

Side A (the distributor) begins by generating n Shamir shares over a Galois Field (GF) with a given prime p. As usual, any t+1 valid shares would enable reconstruction of the secret key in a standard Shamir scheme.

In addition, Side A creates an s × t matrix in GF. Here s represents the maximum number of fake shares and each of the s rows represents a polynomial of degree at most t–1.

For each share, Side A computes a vector:

If the share is valid, the share’s vector is defined as:
r = ⟨ x, p₁(x), p₂(x), …, pₛ(x) ⟩
where x is the share’s x-coordinate and each pᵢ(x) is obtained by evaluating one of the s polynomials.

If the share is fake, the share’s vector is defined as:
r = ⟨ x, r₁, r₂, …, rₛ ⟩
where x is the same x-coordinate, but the subsequent entries r₁, r₂, …, rₛ are chosen at random (thus indistinguishable in distribution from genuine evaluations).

Side A then organizes all of these share vectors as columns of a matrix A.

Finally, Side A transmits both the matrix A and the n Shamir shares to Side B (the receiver).

A key point is that because the x-coordinates are randomly selected, there is no way to tell which columns of matrix A come from valid or fake shares if there are fewer than t+1 valid shares. Also, by deliberately inserting fake shares, Side A makes it impossible for Side B to know the exact number of valid shares unless there are sufficiently many (more than t) valid ones.

On the receiving side, Side B obtains matrix A along with the n Shamir shares. Both parties are aware of the parameters s and t.

Side B then expands the matrix A to form a new matrix, which we call M:

First, Side B extracts the first row of A and calls it the vector v (which contains all the x-coordinates).

Next, Side B creates t new rows by computing v raised to increasing powers. That is, for each i from 0 to t–1, a row is formed by raising every element of v to the power i (denoted v^i). These rows, which form a polynomial basis, are then appended to the top of A. The resulting matrix is named M.

Side B computes the null space of matrix M. If this null space is trivial (only the zero vector), it indicates that there are not enough valid shares for reconstruction.

Using the null space, Side B identifies the indices of the matrix columns where the corresponding entries in the null space basis vectors are nonzero (or, more precisely, are consistently zero at the same positions for fake shares). These indices correspond to the valid shares.

Side B extracts the first elements (x-coordinates) of these identified columns from matrix A.
Side B now knows which share is valid or fake.

Finally, Side B uses these valid shares to reconstruct the secret key K using standard Shamir reconstruction techniques (for example, Lagrange interpolation).

---

### ANALYSIS OF GAUSSIAN ELIMINATION AND PROBABILITY OF FAILURE

The key insight is that the t rows generated from the vector v (i.e., the polynomial basis rows) are used in Gaussian elimination to essentially “cancel out” the contributions of valid shares. In ideal circumstances, if all fake shares were grouped together at one side, the extended matrix would, after elimination, appear as:

C = [ X;  R  0 ]

In this representation, X (a t × m submatrix) corresponds to the t new rows created by raising the x-coordinate vector to successive powers. The zero block, on the other hand, corresponds to the valid share part after elimination. Random fake shares cannot be eliminated into a zero block because they do not originate from the underlying polynomials. Although Side B does not know the secret s polynomials and therefore cannot perform a perfect elimination, it is important to note that Gaussian elimination does not alter the null space of the matrix. This means that the null space computed by Side B will be identical to the null space of the ideal matrix C. With fewer than t+1 valid columns, the receiver cannot determine which columns are valid since the necessary cancellation (that would produce the zero block in the valid share section) does not occur without knowledge of the underlying polynomials.

It follows that the null space of C (and thus of M) contains zeros at the positions corresponding to fake columns. This is why at least t+1 valid (real) columns are necessary in order to obtain a non-trivial null space that Side B can use to reliably filter out the fake shares.

In the borderline case of exactly t valid shares, the top-right t × t submatrix D is invertible. The null space in this case would be trivial. Therefore, having fewer than t+1 valid shares precludes reconstruction—this is an intentional security feature.

A potential source of failure occurs if, by chance, the additional s values in a fake share’s vector exactly match the values that would be produced by the secret s polynomials. In other words, the problem arises when the fake share’s vector accidentally mimics a valid share’s vector. For a single fake share, the probability that all s randomly chosen numbers coincide with those computed by the valid polynomials is (1/p)^s, where p is the prime defining the finite field GF(p).

If there are multiple fake shares, the overall failure probability is given by:
  Failure Probability = 1 - (1 - (1/p)^s)^(number of fake shares)

For instance, if p is 524287 and s = 1, then the chance that a fake share's vector accidentally matches that of a valid share is about 1 in 524287 (roughly 0.000002). For s ≥ 2, the probability becomes virtually negligible. This indicates that, under proper protocol execution, the chance that a fake share’s vector accidentally “passes” as a valid one is extremely low.

It is possible to add an additional verification step—such as reselecting and testing the vector for a fake share—if one wants to further minimize this risk, although this extra measure may slightly compromise the statistical indistinguishability between fake and valid share vectors.

---

### SUMMARY

In summary, the protocol works by embedding a polynomial structure into the share vectors and then using null space analysis to pinpoint valid shares. The Gaussian elimination process in the extended matrix naturally cancels out contributions from valid shares (if enough exist), whereas fake shares, being random, cannot be similarly eliminated. Thus, only when there are at least t+1 valid shares does the null space become non-trivial, allowing for successful reconstruction of the secret key K. Failure under correct protocol execution can occur only when the minimum threshold of t+1 valid shares is not reached or in exceedingly rare cases when a fake share randomly mimics a valid one.
