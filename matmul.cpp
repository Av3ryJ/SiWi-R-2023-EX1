#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
/*
extern "C" {
#include <mkl_cblas.h>
#include <mkl.h>
}
*/
void use_naive(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void use_blas(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void use_Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n, int timing_runs);
void Transposed(double *MatA, double *MatB, double *MatC, int m, int k, int n);
double *transpose(double *mat, int m, int k);
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
        //cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, MatA, m, MatB, n, 0.0, MatC, m );
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
    double *B_transposed = transpose(MatB, m, k);
    naive(MatA, B_transposed, MatC, m, k, n);
}

double *transpose(double *mat, int m, int k) {
    double *mat_t = new double [m*k];
    int blocksize = 16;
    for (; blocksize > 0; blocksize--) {
        if (m%blocksize == 0) break;
    }
    if (blocksize == 0) blocksize = 1;
    int rest = k%blocksize;

    //blocking for majority of matrix
    for (int row = 0; row < m; row+blocksize) { // Achtung! hier muss glaub ich ein '+=' -Maxi
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
            mat_t[((int(k/blocksize)*blocksize) + col)*m + row] = mat[row*k + (int(k/blocksize)*blocksize) + col];
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
    // Matrix A m*k, Matrix B k*n --> Matrix C m*n
    //1. If size >> 512 Divide MatA, MatB and MatC in four sub mats
    std::cout << "Reached Strassen()" << std::endl;
    std::cout << "m, k, n: " << m << k << n << std::endl;
    if (m < 512 || k < 512 || n < 512) {
        // TODO: Do naive or transposed
        return;
    }
    //1.1 Zero pad if necessary

    //1.2 find next 2^n to pad all mats to
    int largest = std::max(m, std::max(k, n));
    int size = 512; // start at 512 bc cannot be smaller
    while (largest > size) size *= 2;
    size /= 2; // size of the four sub-matrices
    // calc needed padding:
    int padA = size-m; // pad MatA unten
    int padBoth = size-m; // pad MatA rechts und B unten
    int padB = size-m; //pad MatB rechts
    std::cout << "Size is: " << size << std::endl;

    //1.2 Divide
    //Only malloc one large array to save time
    double *A11pA22 = new double[size * 14]; // A11 + A22
    double *B11pB22 = A11pA22+(size*1); // B11 + B22
    double *A21pA22 = A11pA22+(size*2); // etc.
    double *B12mB22 = A11pA22+(size*3);
    double *B21mB11 = A11pA22+(size*4);
    double *A11pA12 = A11pA22+(size*5);
    double *A21mA11 = A11pA22+(size*6);
    double *B11pB12 = A11pA22+(size*7);
    double *A12mA22 = A11pA22+(size*8);
    double *B21pB22 = A11pA22+(size*9);
    double *A11 = A11pA22+(size*10);
    double *A22 = A11pA22+(size*11);
    double *B11 = A11pA22+(size*12);
    double *B22 = A11pA22+(size*13);
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

    std::cout << "Matrix and value stuff done" << std::endl;

    for (int row = 0; row < size; ++row){
        off_r = row+size;
        std::cout << "Populating row: " << row << std::endl;
        for (int col = 0; col < size; ++col){
            off_c = col+size;
            //std::cout << col << std::endl;
            //values of mats first to reduce comparisons
            val_A11 = 0; //MatA[col+(row*k)];
            val_A12 = 0; //off_c > k ? 0 : MatA[off_c+(row*k)];
            val_A21 = 0; //off_r > m ? 0 : MatA[col+(off_r*k)];
            val_A22 = 0; //off_c > k || off_r > m ? 0 : MatA[(off_c)+(off_r*k)];
            val_B11 = 0; //MatB[col+(row*n)];
            val_B12 = 0; //off_c > n ? 0 : MatB[off_c+(row*n)];
            val_B21 = 0; //off_r > k ? 0 : MatB[col+(off_r*n)];
            val_B22 = 0; //off_c > n || off_r > k ? MatB[(off_c)+(off_r*n)] : 0;

            A11pA22[col+(row*size)] = val_A11 + val_A22; // Padding only in A12 A21 A22
            B11pB22[col+(row*size)] = val_B11 + val_B22; // Padding only in B12 B21 B22
            A21pA22[col+(row*size)] = val_A21 + val_A22; // offset = size; if (offset+<row/col> larger than <m/k/n> use 0
            B12mB22[col+(row*size)] = val_B12 - val_B22;
            B21mB11[col+(row*size)] = val_B21 - val_B11;
            A11pA12[col+(row*size)] = val_A11 + val_A12;
            A21mA11[col+(row*size)] = val_A21 - val_A11;
            B11pB12[col+(row*size)] = val_B11 + val_B12;
            A12mA22[col+(row*size)] = val_A12 - val_A22;
            B21pB22[col+(row*size)] = val_B21 + val_B22;

            A11[col+(row*size)] = val_A11;
            A22[col+(row*size)] = val_A22;
            B11[col+(row*size)] = val_B11;
            B22[col+(row*size)] = val_B22;
        }
    }

    std::cout << "Populated Mats" << std::endl;

    //Only one new "M"-Array to save time malloc-ing
    double *M1 = new double[size * 7]; // (A11 + A22)(B11 + B22)
    double *M2 = M1+(size * 1); // (A21 + A22)B11
    double *M3 = M1+(size * 2); // A11(B12 - B22)
    double *M4 = M1+(size * 3); // A22(B21 - B11)
    double *M5 = M1+(size * 4); // (A11 + A12)B22
    double *M6 = M1+(size * 5); // (A21 - A11)(B11 + B12)
    double *M7 = M1+(size * 6); // (A12 - A22)(B21 + B22)
    //2. recall self with smaller mats
    Strassen(A11pA22, B11pB22, M1, size, size, size);
    Strassen(A21pA22, B11, M2, size, size, size);
    Strassen(A11, B12mB22, M3, size, size, size);
    Strassen(A22, B21mB11, M4, size, size, size);
    Strassen(A11pA12, B22, M5, size, size, size);
    Strassen(A21mA11, B11pB12, M6, size, size, size);
    Strassen(A12mA22, B21pB22, M7, size, size, size);
    //3. "rebuild" MatC
    for (int row = 0; row < size; ++row){
        off_r = row+size;
        for (int col = 0; col < size; ++col){
            off_c = col+size;
            // C11 = M1 + M4 - M5 + M7
            MatC[col+(row*n)] = M1[col+(row*size)]+M4[col+(row*size)]-M5[col+(row*size)]+M7[col+(row*size)]; // C11 = M1 + M4 - M5 + M7
            if (off_c < n) MatC[off_c+(row*n)] = M3[col+(row*size)]+M5[col+(row*size)]; // C12 = M3 + M5
            if (off_r < m) MatC[col+(off_r*n)] = M2[col+(row*size)]+M4[col+(row*size)]; // C21 = M2 + M4
            if (off_c < n && off_r < m) { // C22 = M1 - M2 + M3 + M6
                MatC[off_c+(off_r*n)] = M1[col+(row*size)]-M2[col+(row*size)]+M3[col+(row*size)]+M6[col+(row*size)];
            }
        }
    }
}