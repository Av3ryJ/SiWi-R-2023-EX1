#include "Timer.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
// Matrix A m*k, Matrix B k*n --> Matrix C m*n

int main(int argc, char* argv[]){
int m;
int k;
int n;

if (argc < 4) {
    std::cout << "Usage: matmul <mat1.in> <mat2.in> <output.out>
}

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

std::string outfile = argv[3];
    
double MatC[m*n];

double time = 100.0;

#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "array" );
#endif

siwir::Timer timer;

   for( int n = 0; n < 1000; ++n )
   {
    timer.reset();
    // the next for is needed bc timer.h is too slow to track small multiplications
    for( int i = 0; i < 10000; ++i ){
        for(int zeilen=0; zeilen<m; ++zeilen){
            for(int spaltenB=0; spaltenB<n; ++spaltenB){
                for(int spaltenA=0; spaltenA<k; ++spaltenA){
                    MatC[zeilen*n+spaltenA]+= MatA[zeilen*k+spaltenA]*MatB[spaltenA*n+spaltenB];
                }
            }
        }
      }
    time = std::min(time, timer.elapsed());
   }

#ifdef USE_LIKWID
   likwid_markerStopRegion( "array" );
#endif

ofstream f (outfile);
f << m << n << std::endl;
for(int i = 0; i<(m*n); ++i) {
    f << MatC[i] << "\n";
}
f.close();
}

/*
matrix Mul:
c11 = a11*b11 + a12*b21
c11
*/
