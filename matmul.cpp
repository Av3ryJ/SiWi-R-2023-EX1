/* in der Konsole:
module load mkl
module load intel
module load openmpi
module load likwid/5.2.2

für Messungen mit Likwid:
likwid-perfctr -C 0 -g FLOPS_DP (oder L2 oder L2CACHE) -m ./matmul matrices/testMatrices/A.in matrices/testMatrices/B.in C.out STD/BLAS/OPT
*/

#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

/*
extern "C" {
#include <mkl_cblas.h>
#include <mkl.h>
}
*/
#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

#define MIN_STRASSEN_SIZE 32
#define TIMING_RUNS 1

void use_naive(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void use_blas(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void use_Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void transpose(double *mat, int m, int k, double *mat_t);
void use_Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n));
void Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n));
void StrassenQuad(double *MatA, double *MatB, double *MatC, int s, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n));
void null_matrix(double *MatC, int size);

int main(int argc, char* argv[]){
    int m;
    int k;
    int n;
    // Matrix A m*k, Matrix B k*n --> Matrix C m*n
    
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
        use_Strassen(MatA, MatB, MatC, m, k, n, Transposed);
    }
    if (var == "STD"){
        use_naive(MatA, MatB, MatC, m, k, n);
    }
    if (var == "BLAS"){
        use_blas(MatA, MatB, MatC, m, k, n);
    }
    if (var == "OPT1"){
        use_Transposed(MatA, MatB, MatC, m, k, n);   //transposed
    }
    if (var == "OPT2"){
        use_Strassen(MatA, MatB, MatC, m, k, n, naive);   //strassen
    }
    if (var == "OPT3"){
        use_Strassen(MatA, MatB, MatC, m, k, n, Transposed);   //strassen + transposed for smaller blocks
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

void use_naive(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    
    double time = 100.0;
    siwir::Timer timer;
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "naiv" );
    #endif

    for( int x = 0; x < TIMING_RUNS; ++x ) {
        timer.reset();
        naive(MatA, MatB, MatC, m, k, n);
        sleep(2);
        time = std::min(time, timer.elapsed());
        if (x != TIMING_RUNS - 1) {
            null_matrix(MatC, m*n);
        }
    }
    time -= 2;
    std::cout << time << " STD" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "naiv" );
         likwid_markerClose();
    #endif
}

void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    // the next for is needed bc timer.h is too slow to track small multiplications
    for(int zeilenC=0; zeilenC<m; ++zeilenC){
        for(int spaltenC=0; spaltenC<n; ++spaltenC){
            for(int spaltenA=0; spaltenA<k; ++spaltenA){ //spaltenA = zeilenB
                MatC[n*zeilenC + spaltenC] += MatA[k*zeilenC + spaltenA] * MatB[n*spaltenA + spaltenC]; // C(i,k) += A(i,j) * B(j,k)
            }
        }
    }
}

void use_blas(double *MatA, double *MatB, double *MatC, int m, int k, int n) {

    double time = 100.0;
    siwir::Timer timer;
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "blas" );
    #endif

    for( int x = 0; x < TIMING_RUNS; ++x ) {
        timer.reset();
        //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, MatA, m, MatB, n, 0.0, MatC, m);
        sleep(2);
        time = std::min(time, timer.elapsed());
        if (x != TIMING_RUNS - 1) {
            null_matrix(MatC, m*n);
        }
    }
    time -= 2;
    std::cout << time << " BLAS" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "blas" );
        likwid_markerClose();
    #endif
}

void use_Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    
    double time = 100.0;
    siwir::Timer timer;
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "transposed" );
    #endif

    for( int x = 0; x < TIMING_RUNS; ++x ) {
        timer.reset();
        Transposed(MatA, MatB, MatC, m, k, n);
        sleep(2);
        time = std::min(time, timer.elapsed());
        if (x != TIMING_RUNS - 1) {
            null_matrix(MatC, m*n);
        }
    }
    time -= 2;
    std::cout << time << " OPT1" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "transposed" );
        likwid_markerClose();
    #endif
}

void Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n) {
    double *B_transposed = new double [k*n];
    transpose(MatB, k, n, B_transposed);
    
    for( int i = 0; i < 1; ++i ){ // TODO: change back to 10000
        for(int zeilenC=0; zeilenC<m; ++zeilenC){
            for(int spaltenC=0; spaltenC<n; ++spaltenC){
                for(int spaltenA=0; spaltenA<k; ++spaltenA){
                    MatC[n*zeilenC + spaltenC] += MatA[k*zeilenC + spaltenA] * B_transposed[k*spaltenC + spaltenA];
                }
            }
        }
    }
    //free B_transposed
    delete [] B_transposed;
}

void transpose(double *mat, int m, int k, double *mat_t) {
    int blocksize = 16;
    for (; blocksize > 0; blocksize--) {
        if (m%blocksize == 0) break;
    }
    if (blocksize == 0) blocksize = 1;
    int rest = k%blocksize;

    //blocking for majority of matrix
    for (int row = 0; row < m; row+=blocksize) {
        for (int col = 0; col < k; col+=blocksize) {
            for (int r = 0; r < blocksize; r++) {
                for (int c = 0; c < blocksize; c++) {
                    mat_t[(col+c)*m + row + r] = mat[(row+r)*k + col + c];
                }
            }
            if ((col + 2*blocksize) > k) col += blocksize;
        }
    }

    //normal for loops for rest of matrix
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < rest; col++) {
            mat_t[((int(k/blocksize)*blocksize) + col)*m + row] = mat[row*k + (int(k/blocksize)*blocksize) + col];
        }
    }
}

void use_Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n)) {
    
    double time = 100.0;
    siwir::Timer timer;
    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "strassen" );
    #endif

    for( int x = 0; x < TIMING_RUNS; ++x ) {
        timer.reset();
        // check if MatA and Matb are 2^n squared matrices
        int size = 2;
        while (m > size) size *=2;
        bool is_squared = (size == m && size == k && size == n);
        // if both are squared we can use StrassenQuad for better performance
        if (is_squared) StrassenQuad(MatA, MatB, MatC, size, function);
        else Strassen(MatA, MatB, MatC, m, k, n, function);
        sleep(2);
        time = std::min(time, timer.elapsed());
        if (x != TIMING_RUNS - 1) {
            null_matrix(MatC, m*n);
        }
    }
    time -= 2;
    std::cout << time << " OPT2/3" << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "strassen" );
        likwid_markerClose();
    #endif
}

void Strassen(double *MatA, double *MatB, double *MatC, int m, int k, int n, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n)) {
    // Matrix A m*k, Matrix B k*n --> Matrix C m*n
    //1. If size >> MIN_STRASSEN_SIZE Divide MatA, MatB and MatC in four sub mats
    
    if (m <= MIN_STRASSEN_SIZE || k <= MIN_STRASSEN_SIZE || n <= MIN_STRASSEN_SIZE) {
        function(MatA, MatB, MatC, m, k, n);
        return;
    }

    //1.2 find next 2^n to pad all mats to
    int largest = std::max(m, std::max(k, n));
    int size = MIN_STRASSEN_SIZE; // start at MIN_STRASSEN_SIZE bc cannot be smaller
    while (largest > size) size *= 2;
    size /= 2; // size of the four sub-matrices
    int sizesquared = size * size;
    //1.2 Divide
    //Only malloc one large array to save time

    double *A11pA22 = new double[sizesquared * 14]; // A11 + A22
    double *B11pB22 = A11pA22+(sizesquared * 1); // B11 + B22
    double *A21pA22 = A11pA22+(sizesquared * 2); // etc.
    double *B12mB22 = A11pA22+(sizesquared * 3);
    double *B21mB11 = A11pA22+(sizesquared * 4);
    double *A11pA12 = A11pA22+(sizesquared * 5);
    double *A21mA11 = A11pA22+(sizesquared * 6);
    double *B11pB12 = A11pA22+(sizesquared * 7);
    double *A12mA22 = A11pA22+(sizesquared * 8);
    double *B21pB22 = A11pA22+(sizesquared * 9);
    double *A11 = A11pA22+(sizesquared * 10);
    double *A22 = A11pA22+(sizesquared * 11);
    double *B11 = A11pA22+(sizesquared * 12);
    double *B22 = A11pA22+(sizesquared * 13);

    //populate with padding if necessary
    int off_r;
    int off_c;
    double val_A11;
    double val_A12;
    double val_A21;
    double val_A22;
    double val_B11;
    double val_B12;
    double val_B21;
    double val_B22;
    bool rgeqm; //row geater or equal m
    bool orgeqm;
    bool rgeqk; //row greater or equal k
    bool orgeqk;
    int index = 0;

    for (int row = 0; row < size; ++row){
        off_r = row+size;
        //compare here and not in iteration
        rgeqm = row >= m;
        orgeqm = off_r >= m;
        rgeqk = row >= k;
        orgeqk = off_r >= k;

        for (int col = 0; col < size; ++col, ++index){
            off_c = col+size;
            //values of mats first to reduce comparisons
            val_A11 = rgeqm  || col >= k   ? 0 : MatA[col+(row*k)];
            val_A12 = rgeqm  || off_c >= k ? 0 : MatA[off_c+(row*k)];
            val_A21 = orgeqm || col >= k   ? 0 : MatA[col+(off_r*k)];
            val_A22 = orgeqm || off_c >= k ? 0 : MatA[(off_c)+(off_r*k)];
            val_B11 = rgeqk  || col >= n   ? 0 : MatB[col+(row*n)];
            val_B12 = rgeqk  || off_c >= n ? 0 : MatB[off_c+(row*n)];
            val_B21 = orgeqk || col >= n   ? 0 : MatB[col+(off_r*n)];
            val_B22 = orgeqk || off_c >= n ? 0 : MatB[(off_c)+(off_r*n)];

            A11pA22[index] = val_A11 + val_A22;
            B11pB22[index] = val_B11 + val_B22;
            A21pA22[index] = val_A21 + val_A22;
            B12mB22[index] = val_B12 - val_B22;
            B21mB11[index] = val_B21 - val_B11;
            A11pA12[index] = val_A11 + val_A12;
            A21mA11[index] = val_A21 - val_A11;
            B11pB12[index] = val_B11 + val_B12;
            A12mA22[index] = val_A12 - val_A22;
            B21pB22[index] = val_B21 + val_B22;

            A11[index] = val_A11;
            A22[index] = val_A22;
            B11[index] = val_B11;
            B22[index] = val_B22;
        }
    }

    //Only one new "M"-Array to save time malloc-ing
    double *M1 = new double[sizesquared * 7]; // (A11 + A22)(B11 + B22)
    double *M2 = M1+(sizesquared * 1); // (A21 + A22)B11
    double *M3 = M1+(sizesquared * 2); // A11(B12 - B22)
    double *M4 = M1+(sizesquared * 3); // A22(B21 - B11)
    double *M5 = M1+(sizesquared * 4); // (A11 + A12)B22
    double *M6 = M1+(sizesquared * 5); // (A21 - A11)(B11 + B12)
    double *M7 = M1+(sizesquared * 6); // (A12 - A22)(B21 + B22)
    //2. recall self with smaller mats
    StrassenQuad(A11pA22, B11pB22, M1, size, function);
    StrassenQuad(A21pA22, B11, M2, size, function);
    StrassenQuad(A11, B12mB22, M3, size, function);
    StrassenQuad(A22, B21mB11, M4, size, function);
    StrassenQuad(A11pA12, B22, M5, size, function);
    StrassenQuad(A21mA11, B11pB12, M6, size, function);
    StrassenQuad(A12mA22, B21pB22, M7, size, function);
    //free A12pA22 + all following matrices
    delete [] A11pA22;
    //3. "rebuild" MatC
    index = 0;
    for (int row = 0; row < size; ++row){
        off_r = row+size;
        for (int col = 0; col < size; ++col, ++index){
            off_c = col+size;
            
            if (row < m && col < n) { // C11 = M1 + M4 - M5 + M7
                MatC[col+(row*n)] = M1[index]+M4[index]-M5[index]+M7[index];
            }
            if (off_c < n) { // C12 = M3 + M5
                MatC[off_c+(row*n)] = M3[index]+M5[index];
            }
            if (off_r < m) { // C21 = M2 + M4
                MatC[col+(off_r*n)] = M2[index]+M4[index];
                if (off_c < n) {
                    MatC[off_c+(off_r*n)] = M1[index]-M2[index]+M3[index]+M6[index];
                }
            }
        }
    }
    //free M1 + all following matrices
    delete [] M1;
}

// This Strassen only works for 2^n square matrices
void StrassenQuad(double *MatA, double *MatB, double *MatC, int s, void (*function)(double *MatA, double *MatB, double *MatC, int m, int k, int n)) {
    // Matrix A m*k, Matrix B k*n --> Matrix C m*n
    //1. If size >> MIN_STRASSEN_SIZE Divide MatA, MatB and MatC in four sub mats
    
    if (s <= MIN_STRASSEN_SIZE) {
        function(MatA, MatB, MatC, s, s, s);
        return;
    }
    //1.1 Zero pad if necessary

    int size = s / 2; // size of the four sub-matrices
    int sizesquared = size * size;

    //1.2 Divide
    //Only malloc one large array to save time

    double *A11pA22 = new double[sizesquared * 14]; // A11 + A22
    double *B11pB22 = A11pA22+(sizesquared *1); // B11 + B22
    double *A21pA22 = A11pA22+(sizesquared*2); // etc.
    double *B12mB22 = A11pA22+(sizesquared*3);
    double *B21mB11 = A11pA22+(sizesquared*4);
    double *A11pA12 = A11pA22+(sizesquared*5);
    double *A21mA11 = A11pA22+(sizesquared*6);
    double *B11pB12 = A11pA22+(sizesquared*7);
    double *A12mA22 = A11pA22+(sizesquared*8);
    double *B21pB22 = A11pA22+(sizesquared*9);
    double *A11 = A11pA22+(sizesquared *10);
    double *A22 = A11pA22+(sizesquared *11);
    double *B11 = A11pA22+(sizesquared *12);
    double *B22 = A11pA22+(sizesquared *13);

    //populate with padding if necessary
    int off_r;
    int off_c;
    double val_A11;
    double val_A12;
    double val_A21;
    double val_A22;
    double val_B11;
    double val_B12;
    double val_B21;
    double val_B22;

    int index = 0;
    for (int row = 0; row < size; ++row){
        off_r = row+size;
        for (int col = 0; col < size; ++col, ++index){
            off_c = col+size;
            //values of mats first to reduce comparisons
            val_A11 = MatA[col+(row*s)];
            val_A12 = MatA[off_c+(row*s)];
            val_A21 = MatA[col+(off_r*s)];
            val_A22 = MatA[(off_c)+(off_r*s)];
            val_B11 = MatB[col+(row*s)];
            val_B12 = MatB[off_c+(row*s)];
            val_B21 = MatB[col+(off_r*s)];
            val_B22 = MatB[(off_c)+(off_r*s)];

            A11pA22[index] = val_A11 + val_A22;
            B11pB22[index] = val_B11 + val_B22;
            A21pA22[index] = val_A21 + val_A22;
            B12mB22[index] = val_B12 - val_B22;
            B21mB11[index] = val_B21 - val_B11;
            A11pA12[index] = val_A11 + val_A12;
            A21mA11[index] = val_A21 - val_A11;
            B11pB12[index] = val_B11 + val_B12;
            A12mA22[index] = val_A12 - val_A22;
            B21pB22[index] = val_B21 + val_B22;

            A11[index] = val_A11;
            A22[index] = val_A22;
            B11[index] = val_B11;
            B22[index] = val_B22;
        }
    }

    //Only one new "M"-Array to save time malloc-ing
    double *M1 = new double[sizesquared * 7]; // (A11 + A22)(B11 + B22)
    double *M2 = M1+(sizesquared * 1); // (A21 + A22)B11
    double *M3 = M1+(sizesquared * 2); // A11(B12 - B22)
    double *M4 = M1+(sizesquared * 3); // A22(B21 - B11)
    double *M5 = M1+(sizesquared * 4); // (A11 + A12)B22
    double *M6 = M1+(sizesquared * 5); // (A21 - A11)(B11 + B12)
    double *M7 = M1+(sizesquared * 6); // (A12 - A22)(B21 + B22)
    //2. recall self with smaller mats
    StrassenQuad(A11pA22, B11pB22, M1, size, function);
    StrassenQuad(A21pA22, B11, M2, size, function);
    StrassenQuad(A11, B12mB22, M3, size, function);
    StrassenQuad(A22, B21mB11, M4, size, function);
    StrassenQuad(A11pA12, B22, M5, size, function);
    StrassenQuad(A21mA11, B11pB12, M6, size, function);
    StrassenQuad(A12mA22, B21pB22, M7, size, function);
    //free A11pA22 + all following matrices
    delete [] A11pA22;
    //3. "rebuild" MatC
    index = 0;
    for (int row = 0; row < size; ++row){
        off_r = row+size;
        for (int col = 0; col < size; ++col, ++index){
            off_c = col+size;
            // C11 = M1 + M4 - M5 + M7
            MatC[col+(row*s)] = M1[index]+M4[index]-M5[index]+M7[index];
            // C12 = M3 + M5
            MatC[off_c+(row*s)] = M3[index]+M5[index];
            // C21 = M2 + M4
            MatC[col+(off_r*s)] = M2[index]+M4[index];
            // C22 = M1 - M2 + M3 + M6
            MatC[off_c+(off_r*s)] = M1[index]-M2[index]+M3[index]+M6[index];
        }
    }
    //free M1 + all following matrices
    delete [] M1;
}
