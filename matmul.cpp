#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
// Matrix A m*k, Matrix B k*n --> Matrix C m*n

extern "C" {
#include <mkl_cblas.h>
#include <mkl.h>
}

int main(int argc, char* argv[]){
int m;
int k;
int n;

//Daten aus Dateien einlesen

std::string dateiA = argv[1];
std::ifstream fileA(dateiA);
fileA>> m;
fileA>>k;
double MatA[m*k]; 
for(int i=0; i<(m*k); i++){
    fileA>>MatA[i];
}
std::string dateiB = argv[2];
std::ifstream fileB(dateiB);
fileB>>k;
fileB>>n;
double MatB[k*n]; 
for(int i=0; i<(k*n); i++){
    fileB>>MatB[i];
}

cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, A, , B, ldb, 0.0, C, ldc );

}