/**
 * Christina Giannoula
 * cgiannoula: christina.giann@gmail.com
 */

#include "./header_files/matrix.h"

#if defined(ONLY_COO)
/**
 * @brief read matrix from input fileName in COO format 
 * @param filename to read matrix (mtx format)
 */
struct COOMatrix *readCOOMatrix(const char* fileName, val_dt s_min, val_dt s_max) { // 1. Same dont touch it (read function of input matrix).

    struct COOMatrix *cooMtx;
    cooMtx = (struct COOMatrix *) malloc(sizeof(struct COOMatrix));
    FILE* fp = fopen(fileName, "r");
    uint32_t rowindx, colindx;
    char *line; 
    char *token; 
    line = (char *) malloc(1000 * sizeof(char));
    int done = false;
    int i = 0;

    while(fgets(line, 1000, fp) != NULL){
        token = strtok(line, " "); //check strtok

        if(token[0] == '%'){
            ;
        } else if (done == false) {
            cooMtx->nrows = atoi(token);
            token = strtok(NULL, " "); 
            cooMtx->ncols = atoi(token);
            token = strtok(NULL, " ");
            cooMtx->nnz = atoi(token);
            printf("[INFO] %s: %u Rows, %u Cols, %u NNZs\n", fileName, cooMtx->nrows, cooMtx->ncols, cooMtx->nnz);
            if((cooMtx->nrows % (8 / byte_dt)) != 0) { // Padding needed
                cooMtx->nrows += ((8 / byte_dt) - (cooMtx->nrows % (8 / byte_dt)));
            }
            if((cooMtx->ncols % (8 / byte_dt)) != 0) {
                cooMtx->ncols += ((8 / byte_dt) - (cooMtx->ncols % (8 / byte_dt))); // Padding needed
            }

            cooMtx->rows = (uint32_t *) calloc((cooMtx->nrows), sizeof(uint32_t));
            cooMtx->nnzs = (struct elem_t *) calloc((cooMtx->nnz+8), sizeof(struct elem_t));
            done = true;
        } else {
            rowindx = atoi(token);
            token = strtok(NULL, " ");
            colindx = atoi(token);
            token = strtok(NULL, " ");

            cooMtx->nnzs[i].rowind = rowindx - 1; // Convert indexes to start at 0
            cooMtx->nnzs[i].colind = colindx - 1; // Convert indexes to start at 0
            cooMtx->nnzs[i].val = (val_dt) (s_min + (s_max - s_min) * drand48());
            i++;
        }
    }

    free(line);
    fclose(fp);
    return cooMtx;
} // 1. End


/**
 * @brief deallocate matrix in COO format 
 * @param matrix in COO format
 */
void freeCOOMatrix(struct COOMatrix *cooMtx) { //3. Don't touch deallocate of input
    free(cooMtx->rows);
    free(cooMtx->nnzs);
    free(cooMtx);
} // 3. End


/** 
 * brief Comparator for Quicksort
 */
int comparator(void *a, void *b) {
    if (((struct elem_t *)a)->rowind < ((struct elem_t *)b)->rowind) {
        return -1;
    } else if (((struct elem_t *)a)->rowind > ((struct elem_t *)b)->rowind) {
        return 1;
    } else {
        return (((struct elem_t *)a)->colind - ((struct elem_t *)b)->colind);
    }
}


/**
 * @brief Sort Input Matrix
 * @param matrix in COO format
 */
void sortCOOMatrix(struct COOMatrix *cooMtx) {

    qsort(cooMtx->nnzs, cooMtx->nnz, sizeof(struct elem_t), comparator);
    
    int prev_row = cooMtx->nnzs[0].rowind;
    int cur_nnz = 0;
    for(unsigned int n = 0; n < cooMtx->nnz; n++) {
        if(cooMtx->nnzs[n].rowind == prev_row)
            cur_nnz++;
        else {
            cooMtx->rows[prev_row] = cur_nnz;
            prev_row = cooMtx->nnzs[n].rowind;
            cur_nnz = 1;
        }
    }
    cooMtx->rows[prev_row] = cur_nnz;

    cur_nnz = 0;
    for(unsigned int r = 0; r < cooMtx->nrows; r++) {
        cur_nnz += cooMtx->rows[r];
    }
    assert(cur_nnz == cooMtx->nnz);
}

#else
/**
 * @brief read matrix from input fileName in COO format 
 * @param filename to read matrix (mtx format)
 */
struct COOMatrix *readCOOMatrix(const char* fileName, val_dt s_min, val_dt s_max) {
    struct COOMatrix *cooMtx;
    cooMtx = (struct COOMatrix *) malloc(sizeof(struct COOMatrix));
    FILE* fp = fopen(fileName, "r");
    uint32_t rowindx, colindx;
    char *line; 
    char *token; 
    line = (char *) malloc(1000 * sizeof(char));
    int done = false;
    int i = 0;

    while(fgets(line, 1000, fp) != NULL){
        token = strtok(line, " ");

        if(token[0] == '%'){
            ;
        } else if (done == false) {
            cooMtx->nrows = atoi(token);
            token = strtok(NULL, " ");
            cooMtx->ncols = atoi(token);
            token = strtok(NULL, " ");
            cooMtx->nnz = atoi(token);
            printf("[INFO] %s: %u Rows, %u Cols, %u NNZs\n", fileName, cooMtx->nrows, cooMtx->ncols, cooMtx->nnz);
            if((cooMtx->nrows % (8 / byte_dt)) != 0) { // Padding needed
                cooMtx->nrows += ((8 / byte_dt) - (cooMtx->nrows % (8 / byte_dt)));
            }
            if((cooMtx->ncols % (8 / byte_dt)) != 0) { // Padding needed
                cooMtx->ncols += ((8 / byte_dt) - (cooMtx->ncols % (8 / byte_dt)));
            }
            cooMtx->rowindx = (uint32_t *) calloc(cooMtx->nnz, sizeof(uint32_t));
            cooMtx->colind = (uint32_t *) calloc(cooMtx->nnz, sizeof(uint32_t));
            cooMtx->values = (val_dt *) calloc(cooMtx->nnz, sizeof(val_dt));
            done = true;
        } else {
            rowindx = atoi(token);
            token = strtok(NULL, " ");
            colindx = atoi(token);
            token = strtok(NULL, " ");

            cooMtx->rowindx[i] = rowindx - 1; // Convert indexes to start at 0
            cooMtx->colind[i] = colindx - 1; // Convert indexes begin at 0
            cooMtx->values[i] = (val_dt) (s_min + (s_max - s_min) * drand48());
            i++;
        }
    }

    free(line);
    fclose(fp);
    return cooMtx;
}


/**
 * @brief deallocate matrix in COO format 
 * @param matrix in COO format
 */
void freeCOOMatrix(struct COOMatrix *cooMtx) {
    free(cooMtx->rowindx);
    free(cooMtx->colind);
    free(cooMtx->values);
    free(cooMtx);
}


/**
 * @brief convert matrix from COO to CSR 
 * @param matrix in COO format
 */
struct CSRMatrix *coo2csr(struct COOMatrix *cooMtx) { // 2. Don't touch

    struct CSRMatrix *csrMtx;
    csrMtx = (struct CSRMatrix *) malloc(sizeof(struct CSRMatrix));

    csrMtx->nrows = cooMtx->nrows;
    csrMtx->ncols = cooMtx->ncols;
    csrMtx->nnz = cooMtx->nnz;
    csrMtx->rowptr = (uint32_t *) calloc((csrMtx->nrows + 2), sizeof(uint32_t));
    csrMtx->colind = (uint32_t *) calloc((csrMtx->nnz + 1), sizeof(uint32_t));
    csrMtx->values = (val_dt *) calloc((csrMtx->nnz + 8), sizeof(val_dt)); // Padding needed

    for(unsigned int i = 0; i < cooMtx->nnz; ++i) {
        uint32_t rowIndx = cooMtx->rowindx[i];
        csrMtx->rowptr[rowIndx]++;
    }

    uint32_t sumBeforeNextRow = 0;
    for(unsigned int rowIndx = 0; rowIndx < csrMtx->nrows; ++rowIndx) {
        uint32_t sumBeforeRow = sumBeforeNextRow;
        sumBeforeNextRow += csrMtx->rowptr[rowIndx];
        csrMtx->rowptr[rowIndx] = sumBeforeRow;
    }
    csrMtx->rowptr[csrMtx->nrows] = sumBeforeNextRow;

    for(unsigned int i = 0; i < cooMtx->nnz; ++i) {
        uint32_t rowIndx = cooMtx->rowindx[i];
        uint32_t nnzIndx = csrMtx->rowptr[rowIndx]++;
        csrMtx->colind[nnzIndx] = cooMtx->colind[i];
        csrMtx->values[nnzIndx] = cooMtx->values[i];
    }

    for(unsigned int rowIndx = csrMtx->nrows - 1; rowIndx > 0; --rowIndx) {
        csrMtx->rowptr[rowIndx] = csrMtx->rowptr[rowIndx - 1];
    }
    csrMtx->rowptr[0] = 0;

    return csrMtx;

} // 2. End
#endif


/**
 * @brief deallocate matrix in CSR format 
 * @param matrix in CSR format
 */
void freeCSRMatrix(struct CSRMatrix *csrMtx) { 
    free(csrMtx->rowptr);
    free(csrMtx->colind);
    free(csrMtx->values);
    free(csrMtx);
}


/**
 * @brief deallocate matrix in BCSR format 
 * @param matrix in BCSR format
 */
void freeBCSRMatrix(struct BCSRMatrix *bcsrMtx) {
    free(bcsrMtx->browptr);
    free(bcsrMtx->bcolind);
    free(bcsrMtx->bval);
    free(bcsrMtx->nnz_per_block);
    free(bcsrMtx);
}


/**
 * @brief deallocate matrix in BCOO format 
 * @param matrix in BCOO format
 */
void freeBCOOMatrix(struct BCOOMatrix *bcooMtx) {
    free(bcooMtx->bind);
    free(bcooMtx->bval);
    free(bcooMtx->nnz_per_block);
    free(bcooMtx);
}


/**
 * @brief print matrix in CSR format
 */
void printCSRMatrix(struct CSRMatrix *A) {

    for(unsigned int rowIndx = 0; rowIndx < A->nrows; ++rowIndx) {
        printf("Row: %d\n", rowIndx);
        for(unsigned int i = A->rowptr[rowIndx]; i < A->rowptr[rowIndx + 1]; ++i) {
            unsigned int colIndx = A->colind[i];
            val_dt value = A->values[i];
            printf("(%d, %f) ", colIndx, value); 
        }
        printf("\n");
    }
}


/**
 * @brief print matrix in BCSR format
 */
void printBCSRMatrix(struct BCSRMatrix *bcsrMtx) {

    printf("Matrix of size %d x %d NNZ %d\n", bcsrMtx->nrows, bcsrMtx->ncols, bcsrMtx->nnz);
    printf("Total Blocks: %ld\n", bcsrMtx->num_blocks);
    printf("Num Block Rows: %d, Num Block Cols: %d of block %ld x %ld\n", bcsrMtx->num_block_rows, bcsrMtx->num_block_cols, bcsrMtx->row_block_size, bcsrMtx->col_block_size);

    for(uint64_t n=0; n<bcsrMtx->num_block_rows; n++) {
        for(uint64_t i=bcsrMtx->browptr[n]; i<bcsrMtx->browptr[n+1]; i++){
            uint64_t j = bcsrMtx->bcolind[i];
            printf("Block (%ld, %ld) NNZ: %d\n", n, j, bcsrMtx->nnz_per_block[i]);
            for(uint64_t r=0; r<bcsrMtx->row_block_size; r++){
                for(uint64_t c=0; c<bcsrMtx->col_block_size; c++) {
                    printf("%f ", bcsrMtx->bval[i * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c]); 
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}


/**
 * @brief print matrix in BCOO format
 */
void printBCOOMatrix(struct BCOOMatrix *bcooMtx) {

    printf("Matrix of size %d x %d NNZ %d\n", bcooMtx->nrows, bcooMtx->ncols, bcooMtx->nnz);
    printf("Total Blocks: %ld\n", bcooMtx->num_blocks);
    printf("Num Block Rows: %d, Num Block Cols: %d of block %ld x %ld\n", bcooMtx->num_block_rows, bcooMtx->num_block_cols, bcooMtx->row_block_size, bcooMtx->col_block_size);

    for(uint64_t n=0; n<bcooMtx->num_blocks; n++) {
        uint32_t i = bcooMtx->bind[n].rowind;
        uint32_t j = bcooMtx->bind[n].colind;
        printf("Block (%d, %d) NNZ: %d\n", i, j, bcooMtx->nnz_per_block[n]);
        for(uint64_t r=0; r<bcooMtx->row_block_size; r++){
            for(uint64_t c=0; c<bcooMtx->col_block_size; c++) {
                printf("%f ", bcooMtx->bval[i * bcooMtx->row_block_size * bcooMtx->col_block_size + r * bcooMtx->col_block_size + c]); 
            }
            printf("\n");
        }
        printf("\n");
    }
}


/**
 * @brief convert matrix from CSR to BCSR format 
 * @param matrix in CSR format
 * Taken from OSKI: A library of automatically tuned sparse matrix kernels
 * http://bebop.cs.berkeley.edu/oski/downloads.html
 */
struct BCSRMatrix *csr2bcsr(struct CSRMatrix *csrMtx, uint32_t row_block_size, uint32_t col_block_size) {
    struct BCSRMatrix *bcsrMtx;
    bcsrMtx = (struct BCSRMatrix *) malloc(sizeof(struct BCSRMatrix));

    uint32_t num_block_rows = (csrMtx->nrows + row_block_size - 1) / row_block_size;
    uint32_t num_block_cols = (csrMtx->ncols + col_block_size - 1) / col_block_size;
    uint32_t num_rows_left = csrMtx->nrows % row_block_size;

    bcsrMtx->nrows = csrMtx->nrows;
    bcsrMtx->ncols = csrMtx->ncols;
    bcsrMtx->nnz = csrMtx->nnz;
    bcsrMtx->row_block_size = row_block_size; //should be divisble by num_rows
    bcsrMtx->col_block_size = col_block_size;
    bcsrMtx->num_block_rows = num_block_rows;
    bcsrMtx->num_block_cols = num_block_cols;
    bcsrMtx->num_rows_left = num_rows_left;

    uint32_t *bAp;
    uint32_t *bAj;
    val_dt *bAx;

    //tmp variables
    uint32_t num_blocks = 0;
    uint32_t *block_count;

    uint32_t *bAp_next;
    uint64_t I, J;
    uint64_t i, j, k, j0, di;
    val_dt a_ij;

    bAp = (uint32_t *) malloc((num_block_rows+2) * sizeof(uint32_t));
    bAp_next = (uint32_t *) malloc((num_block_rows+2) * sizeof(uint32_t));
    block_count = (uint32_t *) malloc((num_block_cols) * sizeof(uint32_t));
    memset(block_count, 0, num_block_cols * sizeof(uint32_t));

    //Phase I: Count the exact number of new blocks to create.
    bAp[0] = 0; 
    if(num_rows_left == 0) {    
        for(I=0; I<num_block_rows; I++) {    
            for(i=I * row_block_size; i < (I+1) * row_block_size; i++) {    
                for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {    
                    j = csrMtx->colind[k];
                    J = j/col_block_size;
                    if(block_count[J] == 0) {    
                        num_blocks++;
                        block_count[J]++;
                    }
                }
            }
            bAp[I+1] = num_blocks;
            for(i = 0; i < num_block_cols; i++)
                block_count[i] = 0;
        }
    } else {
        for(I=0; I<num_block_rows-1; I++) {
            for(i=I * row_block_size; i < (I+1) * row_block_size; i++) {
                for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                    j = csrMtx->colind[k];
                    J = j/col_block_size;
                    if(block_count[J] == 0) {
                        num_blocks++;
                        block_count[J]++;
                    }
                }
            }
            bAp[I+1] = num_blocks;
            for(i = 0; i < num_block_cols; i++)
                block_count[i] = 0;
        }
        for(i = (num_block_rows-1) * row_block_size; i < ((num_block_rows-1) * row_block_size + num_rows_left); i++) {
            for (k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                j = csrMtx->colind[k];
                J = j/col_block_size;
                if (block_count[J] == 0) {
                    num_blocks++;
                    block_count[J]++;
                }
            }
        }
        bAp[num_block_rows] = num_blocks;
        for(i = 0; i < num_block_cols; i++)
            block_count[i] = 0;
    }

    bcsrMtx->num_blocks = num_blocks;
    bAj = (uint32_t *) malloc((num_blocks + 1) * sizeof(uint32_t));
    bAx = (val_dt *) calloc(((num_blocks + 1) * row_block_size * col_block_size), sizeof(val_dt));
    val_dt *blocks = (val_dt *) malloc((row_block_size * col_block_size * num_block_cols) * sizeof(val_dt));
    memset(blocks, 0, row_block_size * col_block_size * num_block_cols * sizeof(val_dt));
    memcpy(bAp_next, bAp, (num_block_rows+1) * sizeof(uint32_t));
    bcsrMtx->nnz_per_block = (uint32_t *) calloc(num_blocks, sizeof(uint32_t));


    //Phase II: Copy all blocks.
    if(num_rows_left == 0) {
        for(I=0; I < num_block_rows; I++) {
            for(i = I * row_block_size, di=0; di<row_block_size; di++, i++) {
                for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                    j = csrMtx->colind[k];
                    J = j / col_block_size;
                    j0 = J * col_block_size;
                    a_ij = csrMtx->values[k];
                    blocks[J * row_block_size * col_block_size + di * col_block_size + j - j0] = a_ij;
                    block_count[J]++;
                }
            }
            for(i = I*row_block_size, di=0; di<row_block_size; di++, i++) {
                for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                    j = csrMtx->colind[k];
                    J = j / col_block_size;
                    j0 = J * col_block_size; 

                    if(block_count[J] > 0) {
                        uint64_t k_next = bAp_next[I]; 
                        bAj[k_next] = J; 
                        memcpy(bAx + k_next * row_block_size * col_block_size, blocks + J * row_block_size * col_block_size, row_block_size * col_block_size * sizeof(val_dt));
                        bAp_next[I]++;
                        assert(bAp_next[I] <= bAp[I+1]);
                        block_count[J] = 0;
                        memset(blocks + J * col_block_size * row_block_size, 0, row_block_size * col_block_size * sizeof(val_dt));
                    }
                }
            }
        }
    } else {
        for(I = 0; I < num_block_rows-1; I++) {
            for(i = I*row_block_size, di=0; di<row_block_size; di++, i++) {
                for(k = csrMtx->rowptr[i]; k<csrMtx->rowptr[i+1]; k++) {
                    j = csrMtx->colind[k];
                    J = j / col_block_size;
                    j0 = J * col_block_size;
                    a_ij = csrMtx->values[k];
                    blocks[J * row_block_size * col_block_size + di * col_block_size + j - j0] = a_ij;
                    block_count[J]++;
                }
            }

            for(i = I*row_block_size, di=0; di < row_block_size; di++, i++) {
                for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                    j = csrMtx->colind[k];
                    J = j / col_block_size;
                    j0 = J * col_block_size;

                    if (block_count[J] > 0) {
                        uint64_t k_next = bAp_next[I];
                        bAj[k_next] = J;
                        memcpy(bAx + k_next * row_block_size * col_block_size, blocks + J * row_block_size * col_block_size, row_block_size * col_block_size * sizeof(val_dt));
                        bAp_next[I]++;
                        assert(bAp_next[I] <= bAp[I+1]);
                        block_count[J] = 0;
                        memset(blocks + J * col_block_size * row_block_size, 0, row_block_size * col_block_size * sizeof(val_dt));
                    }
                }
            }

        }

        for(i = (num_block_rows-1)*row_block_size, di=0; di < num_rows_left; di++, i++) {
            for(k = csrMtx->rowptr[i]; k < csrMtx->rowptr[i+1]; k++) {
                j = csrMtx->colind[k];
                J = j / col_block_size;
                j0 = J * col_block_size;
                a_ij = csrMtx->values[k];
                blocks[J * row_block_size * col_block_size + di * col_block_size + j - j0] = a_ij;
                block_count[J]++;
            }
        }

        for(i = (num_block_rows-1)*row_block_size, di=0; di<num_rows_left; di++, i++) {
            for(k = csrMtx->rowptr[i]; k<csrMtx->rowptr[i+1]; k++) {
                j = csrMtx->colind[k];
                J = j / col_block_size;
                j0 = J * col_block_size;

                if(block_count[J] > 0) {
                    uint64_t k_next = bAp_next[num_block_rows-1];
                    bAj[k_next] = J; 
                    memcpy(bAx + k_next * row_block_size * col_block_size, blocks + J * row_block_size * col_block_size, row_block_size * col_block_size * sizeof(val_dt));
                    bAp_next[num_block_rows-1]++;
                    assert(bAp_next[num_block_rows-1] <= bAp[num_block_rows]);
                    block_count[J] = 0;
                    memset(blocks + J * col_block_size * row_block_size, 0, row_block_size * col_block_size * sizeof(val_dt));
                }
            }
        }

    }


    free(block_count);
    free(blocks);
    free(bAp_next);

    bcsrMtx->browptr = bAp;
    bcsrMtx->bcolind = bAj;
    bcsrMtx->bval = bAx;

    return bcsrMtx;
}


/**
 * @brief convert matrix from BCSR to BCOO format 
 * @param matrix in BCSR format
 */
struct BCOOMatrix *bcsr2bcoo(struct BCSRMatrix *bcsrMtx) {
    struct BCOOMatrix *bcooMtx;
    bcooMtx = (struct BCOOMatrix *) malloc(sizeof(struct BCOOMatrix));

    bcooMtx->nrows = bcsrMtx->nrows;
    bcooMtx->ncols = bcsrMtx->ncols;
    bcooMtx->nnz = bcsrMtx->nnz;
    bcooMtx->num_block_rows = bcsrMtx->num_block_rows;
    bcooMtx->num_block_cols = bcsrMtx->num_block_cols;
    bcooMtx->num_blocks = bcsrMtx->num_blocks;
    bcooMtx->num_rows_left = bcsrMtx->num_rows_left;
    bcooMtx->row_block_size = bcsrMtx->row_block_size;
    bcooMtx->col_block_size = bcsrMtx->col_block_size;

    bcooMtx->bind = (struct bind_t *) malloc(bcooMtx->num_blocks * sizeof(struct bind_t));
    bcooMtx->bval = (val_dt *) malloc((bcooMtx->num_blocks + 1) * bcooMtx->row_block_size * bcooMtx->col_block_size * sizeof(val_dt));
    bcooMtx->nnz_per_block = (uint32_t *) malloc(bcooMtx->num_blocks * sizeof(uint32_t));


    for(uint64_t n=0; n<bcsrMtx->num_block_rows; n++) {
        for(uint64_t i=bcsrMtx->browptr[n]; i<bcsrMtx->browptr[n+1]; i++){
            bcooMtx->bind[i].rowind = n;
            bcooMtx->bind[i].colind = bcsrMtx->bcolind[i];
        }
    }

    memcpy(bcooMtx->bval, bcsrMtx->bval, (bcooMtx->num_blocks + 1) * bcooMtx->row_block_size * bcooMtx->col_block_size * sizeof(val_dt));
    memcpy(bcooMtx->nnz_per_block, bcsrMtx->nnz_per_block, bcooMtx->num_blocks * sizeof(uint32_t));
    return bcooMtx;
}


/**
 * @brief partition function for quickSort 
 */
uint32_t partitionBCSRMatrix(struct BCSRMatrix *bcsrMtx, uint32_t low, uint32_t high) {

    uint32_t pivot = bcsrMtx->bcolind[high]; 
    uint32_t i = low - 1;
    uint32_t temp_ind;
    val_dt temp_val;

    for(uint32_t j = low; j <= high - 1; j++) {
        if(bcsrMtx->bcolind[j] < pivot) {
            i++;
            // swap(i, j)
            temp_ind = bcsrMtx->bcolind[i];
            bcsrMtx->bcolind[i] = bcsrMtx->bcolind[j];
            bcsrMtx->bcolind[j] = temp_ind;
            for(uint32_t r=0; r<bcsrMtx->row_block_size; r++) {
                for(uint32_t c=0; c<bcsrMtx->col_block_size; c++) {
                    temp_val = bcsrMtx->bval[i * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c];
                    bcsrMtx->bval[i * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c] = bcsrMtx->bval[j * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c];
                    bcsrMtx->bval[j * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c] = temp_val;
                }
            }
        }
    }

    // swap(i+1, high)
    temp_ind = bcsrMtx->bcolind[i+1];
    bcsrMtx->bcolind[i+1] = bcsrMtx->bcolind[high];
    bcsrMtx->bcolind[high] = temp_ind;
    for(uint32_t r=0; r<bcsrMtx->row_block_size; r++) {
        for(uint32_t c=0; c<bcsrMtx->col_block_size; c++) {
            temp_val = bcsrMtx->bval[(i+1) * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c];
            bcsrMtx->bval[(i+1) * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c] = bcsrMtx->bval[high * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c];
            bcsrMtx->bval[high * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c] = temp_val;
        }
    }
    return (i+1); 
}

/**
 * @brief quickSort for matrix in BCSR format
 */
void quickSortBCSRMatrix(struct BCSRMatrix *bcsrMtx, int32_t low, int32_t high) {

    if (low < high) {
        uint32_t mid = partitionBCSRMatrix(bcsrMtx, low, high);

        quickSortBCSRMatrix(bcsrMtx, low, mid - 1); 
        quickSortBCSRMatrix(bcsrMtx, mid + 1, high); 
    }
}

/**
 * @brief Sort wrapper for matrix in BCSR format
 */
void sortBCSRMatrix(struct BCSRMatrix *bcsrMtx) {

    for(uint32_t n=0; n<bcsrMtx->num_block_rows; n++) {
        int32_t low = bcsrMtx->browptr[n]; 
        int32_t high = bcsrMtx->browptr[n+1] - 1; 

        quickSortBCSRMatrix(bcsrMtx, low, high);
    }

}

/**
 * @brief count nnz per block in BCSR format matrix
 */
void countNNZperBlockBCSRMatrix(struct BCSRMatrix *bcsrMtx) {
    for(uint32_t n=0; n<bcsrMtx->num_block_rows; n++) {
        for(uint32_t i=bcsrMtx->browptr[n]; i<bcsrMtx->browptr[n+1]; i++){
            for(uint32_t r=0; r<bcsrMtx->row_block_size; r++){
                for(uint32_t c=0; c<bcsrMtx->col_block_size; c++) {
                    if (bcsrMtx->bval[i * bcsrMtx->row_block_size * bcsrMtx->col_block_size + r * bcsrMtx->col_block_size + c] != 0) 
                        bcsrMtx->nnz_per_block[i]++;
                }
            }
        }
    }
}

