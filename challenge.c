#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// =========================================================
// FUNCTION PROTOTYPE
// =========================================================
void sparse__multiply(
    int rows,
    int cols,
    const double* A,
    const double* x,
    int* out_nnz,
    double* values,
    int* col_indices,
    int* row_ptrs,
    double* y
);

// =========================================================
// TODO: USER IMPLEMENTATION
// =========================================================
void sparse_multiply(
    int rows, int cols, const double* A, const double* x,
    int* out_nnz, double* values, int* col_indices, int* row_ptrs,
    double* y
) {
    // Part 1: Extract matrix data into Compressed Sparse Row (CSR) format

    int values_index = 0; // Index for assigning the nonzero elements in the values array
    int nnz = 0; // Counter for each row's nonzero elements
    
    row_ptrs[0] = 0; // Initial index for slicing the rows in A
    
    for (int i = 0; i < rows; i++) // Iterate over the rows of A 
    {
        for (int j = 0; j < cols; j++) // Iterate over the columns of A
        {
            if (fabs(A[i*cols + j]) > 1e-7) // 1e-7 as tolerance for identifying nonzero values
            {
                values[values_index] = A[i*cols + j]; // Save nonzero element
                col_indices[values_index] = j; // Save columns index for the nonzero value
                values_index++; // Increase index counter for values array
                nnz++; // Increase counter for the number of nonzero elements in a row
            }
        }
        row_ptrs[i+1] = nnz; // Save next index for slicing the rows
    }
    *out_nnz = nnz; // Save the final counter for nonzero elements

    // Part 2: Multiply y = A*x, using only the nonzero values in CSR

    // Indexes for creating the slice to identify the current row
    // The multiplication is done only with the nonzero values, no need for multiplying the zero values
    
    // Strategy: slice the rows and use the indices for locating the relevant elements of x for each operation
    // then accumulate the product subresults into the corresponding element of y.
    for (int i = 0; i < rows; i++) // Iterate over the nonzero values
    {
        int index_A = row_ptrs[i]; // Define lower index for slicing the row
        int index_B = row_ptrs[i+1]; // Define upper index for slicing the row

        y[i] = 0.0; // Preliminar value as 0 to avoid garbage before accumulating the answer

        for (int k = index_A; k < index_B; k++) // Iterated between the row slice
        {
            y[i] += values[k] * x[col_indices[k]]; // Accumulate the product into the corresponding element of y
        }
    }
}

// =========================================================
// TEST HARNESS
// =========================================================
int main(void) {
    srand(time(NULL));
    
    const int num_iterations = 100;
    int passed_count = 0;

    for (int iter = 0; iter < num_iterations; ++iter) {
        int rows = rand() % 41 + 5;
        int cols = rand() % 41 + 5;
        double density = 0.05 + (rand() / (double) RAND_MAX) * 0.35;
        
        size_t mat_sz = (size_t) rows * cols;

        double* A = calloc(mat_sz, sizeof(double));
        for (size_t i = 0; i < mat_sz; ++i) {
            if (((double) rand() / RAND_MAX) < density) {
                A[i] = ((double) rand() / RAND_MAX) * 20.0 - 10.0;
            }
        }

        double* values = malloc(mat_sz * sizeof(double));
        int* col_indices = malloc(mat_sz * sizeof(int));
        int* row_ptrs = malloc((rows + 1) * sizeof(int));
        double* x = malloc(cols * sizeof(double));
        double* y_user = malloc(rows * sizeof(double));
        double* y_ref = calloc(rows, sizeof(double));
        int out_nnz = 0;

        for (int i = 0; i < cols; ++i) {
            x[i] = ((double) rand() / RAND_MAX) * 20.0 - 10.0;
        }

        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += A[i * cols + j] * x[j];
            }
            y_ref[i] = sum;
        }

        sparse_multiply(rows, cols, A, x, &out_nnz, values, col_indices, row_ptrs, y_user);

        double max_err = 0.0;
        int passed = 1;
        for (int i = 0; i < rows; ++i) {
            double diff = fabs(y_user[i] - y_ref[i]);
            double tol = 1e-7 + 1e-7 * fabs(y_ref[i]); // Mixed absolute/relative tolerance
            if (diff > tol) {
                max_err = fmax(max_err, diff);
                passed = 0;
            }
        }

        if (passed) {
            passed_count++;
        }

        printf(
            "Iter %2d [%3dx%3d, density=%.2f, nnz=%4d]: %s (Max error: %.2e)\n",
            iter, rows, cols, density, out_nnz, passed ? "PASS" : "FAIL", max_err
        );

        free(A);
        free(values);
        free(col_indices);
        free(row_ptrs);
        free(x);
        free(y_user);
        free(y_ref);
    }

    printf(
        "\n%s (%d/%d iterations passed)\n",
        passed_count == num_iterations ? "All tests passed!" : "Some tests failed.",
        passed_count, num_iterations
    );
           
    return passed_count == num_iterations ? 0 : 1;
}
