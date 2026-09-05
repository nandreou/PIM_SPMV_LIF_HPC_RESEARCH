/**
 * Christina Giannoula
 * cgiannoula: christina.giann@gmail.com
 */


#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "utils.h"

#if defined(ONLY_COO)
/**
 * @brief nnz in COO matrix format
 */
struct elem_t {
    uint32_t rowind;
    uint32_t colind;
    val_dt val;
};

/**
 * @brief COO matrix format
 */
struct COOMatrix {
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nnz;
    uint32_t *rows;
    struct elem_t *nnzs;
};

#else

/**
 * @brief COO matrix format 
 */
struct COOMatrix {
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nnz;
    uint32_t *rowindx;
    uint32_t *colind;
    val_dt *values;
};
#endif

/**
 * @brief CSR matrix format 
 */
struct CSRMatrix {
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nnz;
    uint32_t* rowptr;
    uint32_t* colind;
    val_dt *values;
};

/**
 * @brief BCSR matrix format 
 */
struct BCSRMatrix {
    uint64_t row_block_size;
    uint64_t col_block_size;
    uint32_t num_block_rows;
    uint32_t num_block_cols;
    uint64_t num_blocks;
    uint32_t num_rows_left;
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nnz;
    uint32_t* browptr;  // row pointer
    uint32_t* bcolind;  // column indices
    val_dt* bval;       // nonzeros
    uint32_t* nnz_per_block;  // nnz per block
};

/**
 * @brief BCOO matrix format 
 */
struct BCOOMatrix {
    uint64_t row_block_size;
    uint64_t col_block_size;
    uint32_t num_block_rows;
    uint32_t num_block_cols;
    uint64_t num_blocks;
    uint32_t num_rows_left;
    uint32_t nrows;
    uint32_t ncols;
    uint32_t nnz;
    struct bind_t *bind; // indexes
    val_dt* bval;       // nonzeros
    uint32_t* nnz_per_block;  // nnz per block
};

/**
 * @brief coordinates in BCOO matrix format 
 */
struct bind_t {
    uint32_t rowind;
    uint32_t colind;
};

struct COOMatrix *readCOOMatrix(const char* fileName, val_dt s_min, val_dt s_max);
void freeCOOMatrix(struct COOMatrix *cooMtx);
void freeCSRMatrix(struct CSRMatrix *csrMtx);
void freeBCSRMatrix(struct BCSRMatrix *bcsrMtx);
void freeBCOOMatrix(struct BCOOMatrix *bcooMtx);
void printCSRMatrix(struct CSRMatrix *A);
void printBCSRMatrix(struct BCSRMatrix *bcsrMtx);
void printBCOOMatrix(struct BCOOMatrix *bcooMtx);
struct CSRMatrix *coo2csr(struct COOMatrix *cooMtx);
struct BCSRMatrix *csr2bcsr(struct CSRMatrix *csrMtx, uint32_t row_block_size, uint32_t col_block_size);
struct BCOOMatrix *bcsr2bcoo(struct BCSRMatrix *bcsrMtx);
void sortBCSRMatrix(struct BCSRMatrix *bcsrMtx);
uint32_t partitionBCSRMatrix(struct BCSRMatrix *bcsrMtx, uint32_t low, uint32_t high);
void quickSortBCSRMatrix(struct BCSRMatrix *bcsrMtx, int32_t low, int32_t high);
void countNNZperBlockBCSRMatrix(struct BCSRMatrix *bcsrMtx);
#if defined(ONLY_COO)
int comparator(void *a, void *b);
void sortCOOMatrix(struct COOMatrix *cooMtx);
#endif

#endif

