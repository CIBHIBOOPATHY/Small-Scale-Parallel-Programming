#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
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

	printf("\nMax nnz: %d\n",max_nnz);

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
    int tol=1000;//Set tolerance value
    for(int i=0;i<nrow;i++)
    {
        if(fabs(z[i]-y[i])<tol)
        {
            valid++;
        }
    }
	if(valid==nrow)
    {
        printf("\nThe serial result and parallel result vectors match. Hence solution is Validated\n");
    }
	else
	{
        printf("\nThe serial result and parallel result vectors do not match. Hence the solution is Not validated\n");
    }
}

/*
 * Function: parallel_ELLPACK
 * --------------------------
 *   Multiplies the ELLPACK matrix with a vector by parallelization
 *
 *   *x contains the vector data; *z stores the results obtained during multiplication
 *
 *   returns: None
 *
 */
void parallel_ELLPACK(int *x,double *z)
{
    int i,j;
    double t2;

    #pragma omp parallel for private(i,j) shared(ell,x,z)
    for(int i=0;i<nrow;i++)
    {
        t2=0.0;
        for(int j=0;j<max_nnz;j++)
        {
            int idx=ell.col_idx[i * max_nnz + j];

            if(ell.val[i * max_nnz + j]==0.0)
            {
                continue;
            }
            t2+=ell.val[i * max_nnz + j]*x[idx];
        }
        z[i]=t2;
    }
}

/*
 * Function: ellpack_matrix_vector_multiply
 * ----------------------------------------
 *   Multiplies the ELLPACK  matrix with a vector in both serial and parallel methods
 *
 *   returns: None
 *
 */
void ellpack_matrix_vector_multiply()
{
    // Declare vector to be multiplied with the ellpack matrix
    int *x=(int *)malloc(nrow*sizeof(int));


    // Declare the result vectors of serial and parallel computations
    double *y=(double *)malloc(nrow*sizeof(double)); // Serial result vector
    double *z=(double *)malloc(nrow*sizeof(double)); // Parallel result vector

    // Input data into x vector
    for(int i=0; i<nrow; i++)
    {
        x[i] = i;
    }

    struct timeval start_time, end_time;
    long comp_time=0;
    float avg_comp_time=0.0;

    for(int i=0;i<1000;i++)
    {
        gettimeofday(&start_time, NULL);

        // serial sparse matrix multiplication
        serial_ELLPACK(x,y);

        gettimeofday(&end_time, NULL);

        //Computation time for performing serial matrix vector multiplication is computed
        comp_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - $

        avg_comp_time+=(float)comp_time/1000000.0f;
    }

	avg_comp_time=avg_comp_time/1000.0;

    printf("\nELLLPACK Computational time for serial :%f seconds\n",avg_comp_time);


    double comp_ptime;
    double avg_comp_ptime=0.0;

    for(int i=0;i<1000;i++)
    {

        double start_ptime=omp_get_wtime();

        // parallel sparse matrix multiplication
        parallel_ELLPACK(x,z);

        double end_ptime=omp_get_wtime();


        //Computation time for performing parallel matrix vector multiplication is computed
        comp_ptime = (end_ptime-start_ptime);
        avg_comp_ptime+=comp_ptime;
    }

	avg_comp_ptime=avg_comp_ptime/1000.0;

    printf("\nELLPACK Computational time for parallel :%f seconds\n",avg_comp_ptime);

    //Validation of results
    validate(y,z);
}

int main()
{
    read_matrix("cage4.mtx");

    convertToELLPACK();
    ellpack_matrix_vector_multiply();

    // Free memory
    free(triplets);
    free(ell.row_size);
    free(ell.col_idx);
    free(ell.val);
    return 0;
}



