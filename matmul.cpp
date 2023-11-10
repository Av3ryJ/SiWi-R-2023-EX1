#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

extern "C" {
#include <mkl_cblas.h>
#include <mkl.h>
}

void use_naive(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void use_blas(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void use_Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void use_Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void null_matrix(double *MatC, int size);

int main(int argc, char* argv[]){
    int m;
    int k;
    int n;
    // Matrix A m*k, Matrix B k*n --> Matrix C m*n
    
    // Number of timing runs. Needed, if multiplying small matrices
    int timing_runs = 1; // TODO: change back to 1000
    
    if (argc < 4) {
        std::cout << "Usage: matmul <mat1.in> <mat2.in> <output.out> optional: <VAR>" << std::endl;
        return -1;
    }

    //MatA
    std::ifstream fileA (argv[1]);
    fileA>> m;
    fileA>> k;

    double *MatA = new double[m*k];

    for (int i=0; i<(m*k); i++) {
        fileA>>MatA[i];
    }

    //MatB
    std::ifstream fileB(argv[2]);
    fileB>> k;
    fileB>> n;

    double *MatB = new double[k*n];

    for(int i=0; i<(k*n); i++){
        fileB>>MatB[i];
    }

    //set outfile
    std::string outfile = argv[3];
    //read in VAR if given
    std::string var = "";
    if (argc == 5) {
        var = argv[4];
    }

    //MatC and time
    double *MatC = new double[m*n];
    

    //Decide which implementation to use
    if (var == ""){
        // TODO: Use fastest
    }
    if (var == "STD"){
        use_naive(MatA, MatB, MatC, m, k, n, timing_runs);
    }
    if (var == "BLAS"){
        use_blas(MatA, MatB, MatC, m, k, n, timing_runs);
    }
    if (var == "OPT1"){
        use_Transposed(MatA, MatB, MatC, m, k, n, timing_runs);   //transposed
    }
    if (var == "OPT2"){
        use_Strassen(MatA, MatB, MatC, m, k, n, timing_runs);   //strassen
    }

    std::ofstream fileO (outfile);
    fileO << m << " " << n << std::endl;
    for(int i = 0; i<(m*n); ++i) {
        fileO << MatC[i] << "\n";
    }
    fileO.close();
    
    delete [] MatA;
    delete [] MatB;
    delete [] MatC;
}

void null_matrix(double *MatC, int size){
    for (int i = 0; i < size; ++i) {
        MatC[i] = 0;
    }
}

void use_naive(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs) {
    
    double time = 100.0;
    siwir::Timer timer;
    
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "array" );
    #endif

    for( int x = 0; x < timing_runs; ++x ) {
        timer.reset();
        naive(MatA, MatB, MatC, m, k, n);
        time = std::min(time, timer.elapsed());
        if (x != timing_runs - 1) {
            null_matrix(MatC, m*n);
        }
    }

    std::cout << time << " STD" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "array" );
    #endif
}

void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    // the next for is needed bc timer.h is too slow to track small multiplications
    for( int i = 0; i < 1; ++i ){ // TODO: change back to 10000
        for(int zeilenC=0; zeilenC<m; ++zeilenC){
            for(int spaltenC=0; spaltenC<n; ++spaltenC){
                for(int spaltenA=0; spaltenA<k; ++spaltenA){ //spaltenA = zeilenB
                    MatC[n*zeilenC + spaltenC] += MatA[k*zeilenC + spaltenA] * MatB[n*spaltenA + spaltenC]; // C(i,k) += A(i,j) * B(j,k)
                }
            }
        }
    }
}

void use_blas(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs) {

    double time = 100.0;
    siwir::Timer timer;
    
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "array" );
    #endif

    for( int x = 0; x < timing_runs; ++x ) {
        timer.reset();
        cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, MatA, m, MatB, n, 0.0, MatC, m );
        time = std::min(time, timer.elapsed());
        if (x != timing_runs - 1) {
            null_matrix(MatC, m*n);
        }
    }

    std::cout << time << " STD" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "array" );
    #endif
}

void use_Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs) {
    
    double time = 100.0;
    siwir::Timer timer;
    
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "array" );
    #endif

    for( int x = 0; x < timing_runs; ++x ) {
        timer.reset();
        Transposed(MatA, MatB, MatC, m, k, n);
        time = std::min(time, timer.elapsed());
        if (x != timing_runs - 1) {
            null_matrix(MatC, m*n);
        }
    }

    std::cout << time << " STD" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "array" );
    #endif
}

void Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    double *B_transposed = transpose(MatB);
    naive(MatA, B_transposed, MatC, m, k, n);
}

double *transpose(double *mat, int m, int k) {
    double *mat_t = new double *[m*k];
    int blocksize = 16;
    for (; blocksize > 0; blocksize--) {
        if (m%blocksize == 0) break;
    }
    if (blocksize == 0) blocksize = 1;
    int rest = k%blocksize;

    //blocking for majority of matrix
    for (int row = 0; row < m; row+blocksize) {
        for (int col = 0; col < k; col+blocksize) {
            for (int r = 0; r < blocksize; r++) {
                for (int c = 0; c < blocksize; c++) {
                    mat_t[(col+c)*m + row + r] = mat[(row+r)*k + col + c];
                }
            }
            if ((col + 2*blocksize) > k) col += blocksize;
        }
    }

    //normal for loops for rest of matrix
    for (int row = 0; row < m; row ++) {
        for (int col = 0; col < rest; col++) {
            mat_t[((int(k/blocksize)*blocksize) + col)*m + row] = mat[row*k + (int(k/blocksize)*blocksize) + col]
        }
    }

    return mat_t;
}

void use_Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs) {
    
    double time = 100.0;
    siwir::Timer timer;
    
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "array" );
    #endif

    for( int x = 0; x < timing_runs; ++x ) {
        timer.reset();
        Strassen(MatA, MatB, MatC, m, k, n);
        time = std::min(time, timer.elapsed());
        if (x != timing_runs - 1) {
            null_matrix(MatC, m*n);
        }
    }

    std::cout << time << " STD" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "array" );
    #endif
}

void Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    //1. If size >> 512 Divide MatA, MatB and MatC in four sub mats
        //1.1 Zero pad if neccessary
        //check if padding on rows or cols is neccessary
        //rows?
        if (m%2 != 0) {
            //pad MatA unten
        } 
        if (k%2 != 0) {
            //pad MatA rechts und B unten
        }
        if (n%2 != 0) {
            //pad MatB rechts
        }
        //1.2 Divide
        double *A1122 = 
        double *M1 = Strassen() //(A11 + A22)(B11 + B22)
    //2. recall self with smaller mats
    //3. "rebuild" MatC
}