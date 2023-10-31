#include "Timer.h"

#include <iostream>
#include <vector>

extern "C" {
#include <mkl_cblas.h>
#include <mkl.h>
}

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif


/*
 * useful likwid commands:
 * likwid-topology
 * likwid-perfctr -C 0 -g GROUP -m ./ex01_example
 *   with GROUP = FLOPS_DP, L2CACHE, L2, ... (from likwid-perfctr -a)
 */

int main()
{   
   const int M = 3;
   const int N = 2;
   const int K = 3;
   
   double A[ M * K ] = { 1.0, 2.0, 3.0,
                         2.0, 4.0, 6.0,
                         3.0, 6.0, 9.0 };
   const int lda = 3;
                   
   double B[ K * N ] = { 1.0, 2.0,
                         2.0, 4.0,
                         3.0, 6.0 };
   const int ldb = 2;
                   
   double C[ M * N ] = { -1.0, -1.0,
                         -1.0, -1.0,
                         -1.0, -1.0 };
   const int ldc = 2;

   double time = 100.0;

   mkl_set_num_threads(1);

#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "array" );
#endif
   
   siwir::Timer timer;
   for( int n = 0; n < 1000; ++n )
   {
      timer.reset();
      for( int i = 0; i < 10000; ++i )
         cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, A, lda, B, ldb, 0.0, C, ldc );
      time = std::min(time, timer.elapsed());
   }
#ifdef USE_LIKWID
   likwid_markerStopRegion( "array" );
#endif   
   
    std::cout << "Calculation took " << time << " seconds\nC =\n"
              << " " << C[0] << " " << C[1] << "\n"
              << " " << C[2] << " " << C[3] << "\n"
              << " " << C[4] << " " << C[5] << "\n\n";
         
   ////////////////              
   // C++ vector //     
   ////////////////
              
   time = 100.0;

   std::vector< double > vA{ 1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 3.0, 6.0, 9.0 };
   
   std::vector< double > vB( K * N );
   vB[0] = 1.0; vB[1] = 2.0;
   vB[2] = 2.0; vB[3] = 4.0;
   vB[4] = 3.0; vB[5] = 6.0;
   
   std::vector< double > vC( M * N, -1.0 );
   
#ifdef USE_LIKWID
   likwid_markerStartRegion( "vector" );
#endif
   
   for( int n = 0; n < 1000; ++n )
   {
      timer.reset();
      for( int i = 0; i < 10000; ++i )
         cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, vA.data(), lda, vB.data(), ldb, 0.0, vC.data(), ldc );
      time = std::min(time, timer.elapsed());
   }

#ifdef USE_LIKWID
   likwid_markerStopRegion( "vector" );
   likwid_markerClose();
#endif   
   
    std::cout << "Calculation took " << time << " seconds\nvC =\n"
              << " " << vC[0] << " " << vC[1] << "\n"
              << " " << vC[2] << " " << vC[3] << "\n"
              << " " << vC[4] << " " << vC[5] << "\n\n";   
}
