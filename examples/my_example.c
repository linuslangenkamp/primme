#include <stdlib.h>
#include <stdio.h>

#include "primme_svds.h"   /* header file for PRIMME SVDS too */ 

#ifndef min
#define min(A,B) ((A)<=(B)?(A):(B))
#endif

/**
 * @typedef primme_mvp_fn_t
 * @brief Function pointer for the matrix-vector product (transpose and standard).
 */
typedef void (*primme_mvp_fn_t)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                                int *transpose, primme_svds_params *primme_svds, int *ierr);

/**
 * @enum primme_target_t
 * @brief Specifies the target singular values to find.
 */
typedef enum primme_target_t {
    TOP,   /* Find the largest singular values. */
    LEAST  /* Find the smallest singular values. */
} primme_target_t;

/**
 * @struct primme_result_t
 * @brief Public struct containing the results of the SVD computation.
 * @attention The user is responsible for freeing this struct using `sparse_svd_free`.
 */
typedef struct primme_result_t {
    int rows;          /* Number of rows (m) of the original matrix. */
    int cols;          /* Number of columns (n) of the original matrix. */
    int target_size;   /* Number of singular values found. */
    double *svals;     /* Array with the computed singular values. */
    double *svecs;     /* Array with the computed singular vectors. */
    double *rnorms;    /* Array with the computed residual norms. */
} primme_result_t;

/**
 * @struct primme_internal_info_t
 * @brief Encapsulates SVD problem data and results.
 * Contains matrix dimensions, config, and result pointers for SVD computations.
 * @note This struct is a private implementation detail, managed via `void*` handle.
 */
typedef struct primme_internal_info_t {
    primme_result_t result;         /* The results of the SVD. */
    primme_svds_params primme_svds; /* PRIMME's internal state. */
} primme_internal_info_t;

/**
 * @brief Allocates and initializes an SVD computation handle.
 * @param rows [in] The number of rows of the matrix.
 * @param cols [in] The number of columns of the matrix.
 * @param target_size [in] The number of singular values to compute.
 * @param linear_operator [in] The callback function for the matrix-vector product.
 * @return A `void*` handle to the internal SVD state, or `NULL` if initialization fails.
 */
static void* sparse_svd_allocate(int rows, int cols, int target_size, primme_mvp_fn_t linear_operator) {
    primme_internal_info_t *info = (primme_internal_info_t*)malloc(sizeof(primme_internal_info_t));
    primme_svds_initialize(&info->primme_svds);
    info->primme_svds.m = rows;
    info->primme_svds.n = cols;
    info->primme_svds.numSvals = target_size;
    info->primme_svds.matrixMatvec = linear_operator;

    info->result.rows = rows;
    info->result.cols = cols;
    info->result.target_size = target_size;
    info->result.svals = (double *) malloc(info->primme_svds.numSvals * sizeof(double));
    info->result.svecs = (double *) malloc((info->primme_svds.n + info->primme_svds.m) * info->primme_svds.numSvals * sizeof(double));
    info->result.rnorms = (double *) malloc(info->primme_svds.numSvals * sizeof(double));

    if (!info->result.svals || !info->result.svecs || !info->result.rnorms) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(info->result.svals);
        free(info->result.svecs);
        free(info->result.rnorms);
        return NULL;
    }

    return info;
}

/**
 * @brief Deallocates all memory associated with the SVD computation handle.
 * @param handle [in] The opaque handle returned by `sparse_svd_allocate`.
 */
static void sparse_svd_free(void* handle) {
    primme_internal_info_t *info = (primme_internal_info_t *)handle;
    primme_svds_free(&info->primme_svds);
    free(info->result.svals);
    free(info->result.svecs);
    free(info->result.rnorms);
    free(info);
}

/**
 * @brief Prints the SVD results to standard output.
 * @param results [in] A pointer to a `primme_result_t` struct containing the SVD results.
 */
static void print_svd_results(primme_result_t* results) {
    if (!results) {
        fprintf(stderr, "Error: primme_result_t pointer is NULL.\n");
        return;
    }

    /* Reporting */
    fprintf(stdout, "\n--- SVD Results ---\n");
    for (int i = 0; i < results->target_size; i++) {
        fprintf(stdout, "sigma[%d]: %-22.15E rnorm: %-22.15E\n", i + 1,
            results->svals[i], results->rnorms[i]);
    }
    fprintf(stdout, "All target_size = %d singular triplets found\n", results->target_size);

    /* Print the largest singular vector (first triplet) */
    fprintf(stdout, "\nLargest left singular vector:\n");
    for (int i = 0; i < results->rows; i++) {
        fprintf(stdout, "%-22.15E\n", results->svecs[i]);
    }

    fprintf(stdout, "\nLargest right singular vector:\n");
    for (int i = 0; i < results->cols; i++) {
        fprintf(stdout, "%-22.15E\n", results->svecs[results->rows + i]);
    }
    fprintf(stdout, "---------------------\n");
}

/**
 * @brief Performs the singular value decomposition.
 * @param handle [in] The opaque handle returned by `sparse_svd_allocate`.
 * @param target [in] Specifies whether to find the `TOP` or `LEAST` singular values.
 * @return A pointer to a `primme_result_t` struct on success, or `NULL` on error.
 */
static primme_result_t* sparse_svd(void* handle, primme_target_t target) {
    primme_internal_info_t *info = (primme_internal_info_t *)handle;

    /* some default values for now */
    double eps = 1e-14;

    /* Set problem matrix */
    info->primme_svds.matrix = NULL; /* user_data */
    info->primme_svds.eps = eps;     /* ||r|| <= eps * ||matrix|| */
    info->primme_svds.target = (target == TOP ? primme_svds_largest : primme_svds_smallest);

    primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
                           PRIMME_DEFAULT_METHOD, &info->primme_svds);

    info->primme_svds.printLevel = 3;

    /* Display PRIMME SVDS configuration struct (optional) */
    primme_svds_display_params(info->primme_svds);

    /* Call primme_svds  */
    int ret = dprimme_svds(info->result.svals, info->result.svecs, info->result.rnorms, &info->primme_svds);

    if (ret != 0) {
        /* error case */
        fprintf(stderr, "Error: primme_svds returned with nonzero exit status: %d \n",ret);
        return NULL;
   }

    return &info->result;
}

static void LinearOperator(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds, int *err);

int main (int argc, char *argv[]) {
    void *handle = sparse_svd_allocate(20000, 20000, 10, LinearOperator);
    primme_result_t* res = sparse_svd(handle, LEAST);
    //print_svd_results(res);
    sparse_svd_free(handle);
}

/* lauchli block matrix-vector product, y = a * x (or y = a^t * x), where

   - x, input dense matrix of size primme_svds.n (or primme_svds.m) x blocksize;
   - y, output dense matrix of size primme_svds.m (or primme_svds.n) x blocksize;
   - a, lauchli matrix of dimensions primme_svds.m x (primme_svds.m+1) with this form:

        [ 1  1  1  1  1 ...   1 ],  ei = 1 - (1 - mu)*i/(min(m,n) - 1)
        [e0  0  0  0  0 ...   0 ]
        [ 0 e1  0  0  0 ...   0 ]
         ...
        [ 0  0  0  0  0 ... en-1]
*/

/**
 * @brief Implements the matrix-vector products for the given matrix.
 * It operates on blocks of vectors for improved performance and 
 * in general looks like this (depending on input):
 *                             Y := A * X,   for transpose = 0
 *                          or Y := A^T * X, for transpose = 1
 *
 * @attention get column i of x: (double *)x + (*ldx) * i;
 * @attention get column i of y: (double *)y + (*ldy) * i;
 *
 * @param x [in] Input dense matrix of vectors `X`.
 * @param ldx [in] Leading dimension of the input matrix `X`.
 * @param y [out] Output dense matrix of vectors `Y`.
 * @param ldy [in] Leading dimension of the output matrix `Y`.
 * @param blockSize [in] Number of vectors in the current block, number of columns of the X matrix.
 * @param transpose [in] Flag indicating if the transpose is applied (0 for A*x, 1 for A^T*x).
 * @param primme_svds [in] PRIMME configuration struct.
 * @param err [out] Error status; must be set to 0 on success.
 */
void LinearOperator(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                    int *transpose, primme_svds_params *primme_svds, int *err) {
   int i;            /* vector index, from 0 to *blockSize-1 */
   int j;
   int min_m_n = min(primme_svds->m, primme_svds->n);
   double *xvec;     /* pointer to i-th input vector x */
   double *yvec;     /* pointer to i-th output vector y */
   void* matrix = (void*)primme_svds->matrix;

   double mu = 1e-1;

   if (*transpose == 0) { /* Do y <- A * x */
      for (i=0; i<*blockSize; i++) {
         xvec = (double *)x + (*ldx)*i;
         yvec = (double *)y + (*ldy)*i;
         yvec[0] = 0;
         for (j=0; j<primme_svds->n; j++) {
            yvec[0] += xvec[j];
         }
         for (j=1; j<primme_svds->m; j++) {
            yvec[j] = j-1<primme_svds->n ? xvec[j-1]*(1.0 - (1.0 - mu)*(j-1)/(min_m_n - 1)) : 0.0;
         }
      }
   } else { /* Do y <- A^t * x */
      for (i=0; i<*blockSize; i++) {
         xvec = (double *)x + (*ldx)*i;
         yvec = (double *)y + (*ldy)*i;
         for (j=0; j<primme_svds->n; j++) {
            yvec[j] = xvec[0];
            if (j+1 < primme_svds->m) yvec[j] += xvec[j+1]*(1.0 - (1.0 - mu)*j/(min_m_n - 1));
         }
      }
   }
   *err = 0;
}
