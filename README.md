# Small-Scale-Parallel-Programming
This project develops a sparse matrix-vector product kernel in C (CSR,ELLPACK formats), parallelized using OpenMP and CUDA for computing y←Ax on Crescent at Cranfield University. Accuracy is verified against a serial implementation. Performance can be analyzed and compared between OpenMP and CUDA versions.

The test matrices are sourced from the Suite Sparse Matrix Collection from the website https://sparse.tamu.edu/.
