# sparse_matrix_csr_challenge
Final submission for a challenge provided to implement a function to execute Compressed Sparse Row (CSR) in C for sparse matrices.

The solution implemented for the function is divided in two parts, according to the instructions:
- Part 1: Compressing the matrix into CSR, by saving the nonzero values, their column indices and the number of nonzero elements per row.
- Part 2: Implementing the multiplication y = A*x, with A being a matrix and x a vector, using only the CSR arrays (with no unnecessary multiplications by zero elements).
