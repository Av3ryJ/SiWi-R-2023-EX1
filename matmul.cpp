#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
// Matrix A m*k, Matrix B k*n --> Matrix C m*n

void naive(double *MatA, double *MatB, double *MatC, int m, int k, int n);
void null_matrix(double *MatC, int size);

int main(int argc, char* argv[]){
    int m;
    int k;
    int n;

    if (argc < 4) {
        std::cout << "Usage: matmul <mat1.in> <mat2.in> <output.out>" << std::endl;
        return -1;
    }

    //MatA
    std::string temp;
    //Daten aus Dateien einlesen
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
    //MatC and time
    double *MatC = new double[m*n];
    double time = 100.0;
    siwir::Timer timer;

    #ifdef USE_LIKWID
        likwid_markerInit();
        likwid_markerStartRegion( "array" );
    #endif


    //naive implementation
    int limit = 2;              // TODO: change back to 1000
    for( int x = 0; x < limit; ++x ) {
        timer.reset();
        naive(MatA, MatB, MatC, m, k, n);
        time = std::min(time, timer.elapsed());
        if (x != limit - 1) {
            null_matrix(MatC, m*n);
        }
    }

    std::cout << "naive algorithm time: " << time << std::endl;

    #ifdef USE_LIKWID
        likwid_markerStopRegion( "array" );
    #endif

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
