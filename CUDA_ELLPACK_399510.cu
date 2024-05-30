#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include<string.h>
#include<float.h>
#include<math.h>

// Structure to represent the matrix input in triplet format
typedef struct {
    int row;
    int col;
    double val;
} Triplet;

// Structure to represent an ELLPACK matrix
typedef struct {
    int *row_size;
    int *col_idx;
    double *val;

} ELLPACKMatrix;

//Matrix information are initialised
int nrow,ncol,nnz,max_nnz;

//Structure creation for matrix data in triplet and ELLPACK formats
Triplet *triplets;
ELLPACKMatrix ell;


/*
 * Function: compare
 * -----------------
 *   Compares row indices of triplet data
 *
 *   returns: 0 if row indices of triplets x and y match;
 *            1 if row index of triplet x greater than row index of y
 *           -1 if row index of triplet x lesser than row index of y
 */
int compare(const void *x, const void *y)
{
    Triplet *tripletsx=(Triplet*)x;
    Triplet *tripletsy=(Triplet*)y;

    if(tripletsx->row > tripletsy->row)
    {
        return 1;
    }
	else if(tripletsx->row < tripletsy->row)
    {
        return -1;
    }
	else
	{
        return 0;
    }
}

/*
 * Function: read_matrix
 * ---------------------
 *   Reads the matrix market file and processes the type of matrix
 *   and stores the data in a triplet structure
 *   
 *   filename: Name of matrix market file to be read
 *   
 *   returns: None
 */
 void read_matrix(char* filename)
{
    char line[512];
    FILE* file=fopen(filename,"r");

    // Read the first line to identify the type of matrix
    fgets(line,sizeof(line),file);

    char *token;
    const char delimiter[2]=" ";

    //get the fourth and fifth word
    token = strtok(line,delimiter);

	int k=0;
    
    int mtype=0; // if 0 - real ; 1 - pattern;
    int mformat=0; // if 0 - general; if 1 - symmetric;
    while(token != NULL)
    {
        k++;
        if(k == 4)
        {
            if(strncmp(token, "pattern", 7) == 0){
                mtype=1;
            }
        } 
        else if(k == 5) 
        {
            if(strncmp(token, "symmetric", 9) == 0){
                mformat=1;
            }

        }
        token = strtok(NULL, delimiter);
    }
	//Read comments
    do{
        fgets(line,sizeof(line),file);
    }while(line[0]=='%');

    // Read basic matrix information
	sscanf(line,"%d %d %d",&nrow,&ncol,&nnz);

    //memory allocation for triplet storage
    if(mformat==1)
    {
        triplets=(Triplet*)malloc(nnz*2*sizeof(Triplet));
	}
    else
	{
        triplets=(Triplet*)malloc(nnz*sizeof(Triplet));
    }


    //Pattern matrix check and read the data
    if(mtype==1)
    {
        for(int i=0;i<nnz;i++)
        {
            fgets(line,sizeof(line),file);
            sscanf(line,"%d %d",&triplets[i].row,&triplets[i].col);
            triplets[i].val=1.0; //store 1.0 for all (row,col) values
            triplets[i].row--; // correction of row index from matrix file
            triplets[i].col--; // correction of colummn index from matrix file
        }
	}
    else
    {
        for(int i=0;i<nnz;i++)
        {
            fgets(line,sizeof(line),file);
            sscanf(line,"%d %d %lf",&triplets[i].row,&triplets[i].col,&triplets[i].val);
            triplets[i].row--;
            triplets[i].col--;
        }
    }

    int index=nnz;

    //Symmetric matrix check to reconstruct the full symmetric matrix
    if(mformat==1)
    {
        for(int i=0;i<nnz;i++)
        {
            if(triplets[i].row!=triplets[i].col)
            {
                triplets[index].row=triplets[i].col;
                triplets[index].col=triplets[i].row;
                triplets[index].val=triplets[i].val;
                index+=1;
            }
        }
    }

	nnz=index;
    qsort(triplets,nnz,sizeof(Triplet),compare);
    fclose(file);
}

/*
 * Function: convertToELLPACK
 * --------------------------
 *   Converts the triplet matrix structure to ELLPACK format
 *
 *   returns: None
 *
 */
 void convertToELLPACK()
{
    //Calculation of maximum number of non-zeroe values in a row
    max_nnz = 0;
    int* nnz_count =(int *)calloc(nrow, sizeof(int));

    for (int i = 0; i < nnz; i++)
    {
        nnz_count[triplets[i].row]++;
        max_nnz = (max_nnz > nnz_count[triplets[i].row]?max_nnz:nnz_count[triplets[i].row]);
    }


    // Allocate memory for ELLPACK format

    ell.col_idx = (int*) malloc(nrow * max_nnz * sizeof(int));
    ell.val = (double*) malloc(nrow * max_nnz * sizeof(double));
    ell.row_size = (int*) calloc(nrow,sizeof(int));

    printf("\nMemory allocated\n");

    // Convert triplets to ELLPACK format
    for (int i = 0; i < nnz; i++)
    {
        int row = triplets[i].row;
        int col = triplets[i].col;
        int idx = ell.row_size[row];

        ell.val[row * max_nnz + idx] = triplets[i].val;
        ell.col_idx[row * max_nnz + idx] = col;
        ell.row_size[row]++;
    }

	// Free memory
    free(nnz_count);
}

/*
 * Function: serial_ELLPACK
 * ------------------------
 *   Multiplies the triplet matrix with a vector in serial
 *
 *   *x contains the vector data; *y stores the results obtained during multiplication
 *   returns: None
 *
 */
void serial_ELLPACK(int *x,double *y)
{
    for(int i=0;i<nrow;i++)
    {
        double t=0.0;
        for(int j=0;j<max_nnz;j++)
        {
            int idx=ell.col_idx[i * max_nnz + j];

            if(ell.val[i * max_nnz + j]==0.0)
            {
                continue;
            }
            t+=ell.val[i * max_nnz + j]*x[idx];
        }
        y[i]=t;
    }
}

/*
 * Function: validate
 * ------------------
 *   Validation of parallel multiplication result vector by comparing it with serial
 *   multiplication result vector
 *
 *   *y contains the serial results; *z contains the parallel results
 *
 *   returns: None
 *
 */
void validate(double *y,double *z)
{
    int valid=0; // Stores the number of matched values
    int tol=1000;
    for(int i=0;i<nrow;i++)
    {
        if(fabs(z[i]-y[i])<tol)
        {
            valid++;
        }
    }
	if(valid==nrow)
    {
        printf("\nValidated\n");
    }
	else
	{
        printf("\nNot validated\n");
    }
}

/*
 * Function: parallel_ELLPACK
 * --------------------------
 *   Parallelization of serial_ELLPACK function
 *
 *   returns: None
 *
 */
 __global__ void parallel_ELLPACK(int* col_idx, double* val, int max_nnz, int nrow, int* x, double* z)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < nrow)
    {
        double t = 0.0;
        for(int j=0; j<max_nnz; j++)
        {
            int idx = col_idx[i * max_nnz + j];
            if(val[i * max_nnz + j] == 0.0)
            {
                continue;
            }
            t += val[i * max_nnz + j] * x[idx];
        }
        z[i] = t;
    }
}

/*
 * Function: ellpack_matrix_vector_multiply
 * -----------------------------------------
 *   Calls the kernels for serial and parallel matrix vector dot product
 *
 *   returns: None
 *
 */
    // initialise vector to be multiplied with the ellpack matrix
    int* x = (int*)malloc(nrow * sizeof(int));


    // initialise the result vectors of serial and parallel computations
    double* y = (double*)malloc(nrow * sizeof(double)); // Serial result vector
    double* z = (double*)malloc(nrow * sizeof(double)); // Parallel result vector

    // Input data into x vector
    for(int i=0; i<nrow; i++)
    {
        x[i] = i;
    }

	// allocate memory on the device
    int* d_col_idx;
    double* d_val;
    int* d_x;
    double* d_z;

    cudaMalloc((void**)&d_col_idx, nrow * max_nnz * sizeof(int));
    cudaMalloc((void**)&d_val, nrow * max_nnz * sizeof(double));
    cudaMalloc((void**)&d_x, nrow * sizeof(int));
    cudaMalloc((void**)&d_z, nrow * sizeof(double));

    // copy matrix data from host to device
    cudaMemcpy(d_col_idx, ell.col_idx, nrow * max_nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, ell.val, nrow * max_nnz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, nrow * sizeof(int), cudaMemcpyHostToDevice);

    // kernel parameters
    int block_size = 256;
    int num_blocks = (nrow + block_size - 1) / block_size;

    float avg_ptime = 0;
    //launch kernel
    for(int i=0;i<1000;i++)
    {
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
        parallel_ELLPACK<<<num_blocks, block_size>>>(d_col_idx, d_val, max_nnz, nrow, d_x, d_z);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float ptime = 0;
        cudaEventElapsedTime(&ptime, start, stop);
        float seconds = ptime / 1000.0f;
        avg_ptime+=seconds;
    }
	printf("\nELLPACK Computational Time for Parallel Kernel:%f seconds\n",avg_ptime/1000);

    // copy results from device to host
    cudaMemcpy(z, d_z, nrow * sizeof(double), cudaMemcpyDeviceToHost);

    // serial sparse matrix multiplication
    serial_ELLPACK(x, y);

    // Validation of results
    validate(y, z);

    // free memory on the device
    cudaFree(d_col_idx);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_z);

    // free memory on the host
    free(x);
    free(y);
    free(z);
}

int main()
{
    //Read matrix
    read_matrix("/scratch/s399510/cage4.mtx");

    convertToELLPACK();
    ellpack_matrix_vector_multiply();

    // Free memory
    free(triplets);
    free(ell.row_size);
    free(ell.col_idx);
    free(ell.val);
    return 0;
}







