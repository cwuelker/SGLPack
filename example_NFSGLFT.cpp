//
//        SGLPack 1.1 - Fast spherical Gauss-Laguerre Fourier transforms        
//   
//  
//   Contact: Christian Wuelker
//            christian.wuelker@gmail.com
//  
//   Copyright 2019  Christian Wuelker, Juergen Prestin
//
//
//     SGLPack is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 3 of the License, or
//     (at your option) any later version.
//  
//     SGLPack is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//     GNU General Public License for more details.
//  
//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//  
//  
//   Commercial use is absolutely prohibited.
//  
//   See the accompanying LICENSE file for details.

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>

#include <nfsglft.h>

int main ()
{
   srand (time (NULL));
   
   const int B = 4;
   const int M = 40000000;
   
   double r_max = 5.0;
   
   const int nfft_sigma = 2;
   const int nfft_q = 5;
   
   const int n_runs = 1;
   
// ----------------------------------------------------------------------------------------------------    
   
   printf ("\n>> Test: NFSGLFT\n");     
   printf (">> r_max: %f...\n", r_max);
   
   //std::ofstream abs_mean ("./abs_mean.txt", std::ofstream::out);
   //std::ofstream abs_stdv ("./abs_stdv.txt", std::ofstream::out);   
   
   //std::ofstream rel_mean ("./rel_mean.txt", std::ofstream::out);
   //std::ofstream rel_stdv ("./rel_stdv.txt", std::ofstream::out);
   
   //abs_mean << std::scientific;
   //abs_stdv << std::scientific;
   
   //rel_mean << std::scientific;
   //rel_stdv << std::scientific;
   
   for ( r_max = 5.0; r_max <= 5.0; r_max += 0.5 )   
   {
      
// ----------------------------------------------------------------------------------------------------      
   
   std::vector<double> runtime_fast_vec;
   std::vector<double> runtime_slow_vec;
   
   std::vector<double> max_abs_error_vec;
   std::vector<double> max_rel_error_vec;
   
   //std::vector<double> max_abs_weight_error_vec;
   //std::vector<double> max_rel_weight_error_vec;   
   
   for ( int run = 0; run < n_runs; ++run )
   {
      
   printf (">> run %d...\n", run);   
   
// ----------------------------------------------------------------------------------------------------
   
   const double pi = 3.141592653589793238462643383279502884197169399;
   
   double*** f_hat_real = new double**[B];
   double*** f_hat_imag = new double**[B];
   
   for ( int n = 1; n <= B; ++n )
   {
      f_hat_real[n - 1] = new double*[n];
      f_hat_imag[n - 1] = new double*[n];
   
      for ( int l = 0; l < n; ++l )
      {
         f_hat_real[n - 1][l] = new double[2 * l + 1];
         f_hat_imag[n - 1][l] = new double[2 * l + 1];
      }
   }
   
   //printf (">> generating (random) SGL Fourier coefficients... ");
   
   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         for ( int m = - l; m <= l; ++m )
         {              
            f_hat_real[n - 1][l][m + l] = (double)(((rand () % 200000) - 100000) / 100000.0);// * exp (- 0.5 * (2 * n - l - 2));
            f_hat_imag[n - 1][l][m + l] = (double)(((rand () % 200000) - 100000) / 100000.0);// * exp (- 0.5 * (2 * n - l - 2));
         }
      }
   }
   
   //printf ("done.\n");
   
   //printf (">> generating %d random evaluation nodes... ", M);

   double** x = new double*[M];
   double* r = new double[M];

   for( int i = 0; i < M; ++i )
   {
      x[i] = new double[3];
      
      long double xi_0, xi_1, xi_2;
      
      do
      {
         xi_0 = (long double)(2.0 * r_max * ((rand () % 100000 - 50000) / 100000.0));
         xi_1 = (long double)(2.0 * r_max * ((rand () % 100000 - 50000) / 100000.0));
         xi_2 = (long double)(2.0 * r_max * ((rand () % 100000 - 50000) / 100000.0));
         
         x[i][0] = sqrt ((double)(xi_0 * xi_0 + xi_1 * xi_1 + xi_2 * xi_2));
         x[i][1] = acos ((double)(xi_2 / x[i][0]));
         x[i][2] = atan2 ((double)xi_1, (double)xi_0);
         
         r[i] = x[i][0];
         
      } while ( r[i] > r_max );
   }
   
   //printf ("done.\n");
   
   double* f_real_ground_truth = new double[M];
   double* f_imag_ground_truth = new double[M];
   
   double* f_real = new double[M];
   double* f_imag = new double[M];
   
// ----------------------------------------------------------------------------------------------------   
   
   //printf (">> running NDSGLFT... ");
   
// ----------------------------------------------------------------------------------------------------

   clock_t begin = clock ();
   
// ----------------------------------------------------------------------------------------------------
   
   ndsglft (B, M, x, f_real_ground_truth, f_imag_ground_truth, f_hat_real, f_hat_imag);
   
// ----------------------------------------------------------------------------------------------------   
   
   clock_t end = clock ();
   
// ----------------------------------------------------------------------------------------------------

   //printf ("done.\n");
   
   double elapsed_msecs = 1000.0 * (double)(end - begin) / CLOCKS_PER_SEC;
   
   if ( run >= 0) runtime_slow_vec.push_back (elapsed_msecs);   
   
   //printf (">> computation time: %f ms\n", elapsed_msecs);
   
// ----------------------------------------------------------------------------------------------------     
   
   //printf (">> running NFSGLFT (sigma: %d, q = %d)... ", nfft_sigma, nfft_q);
   
// ----------------------------------------------------------------------------------------------------

   begin = clock ();
   
// ----------------------------------------------------------------------------------------------------      
   
   nfsglft (B, M, x, f_real, f_imag, f_hat_real, f_hat_imag, nfft_sigma, nfft_q);
   
// ----------------------------------------------------------------------------------------------------   
   
   end = clock ();
   
// ----------------------------------------------------------------------------------------------------

   //printf ("done.\n");
   
   elapsed_msecs = 1000.0 * (double)(end - begin) / CLOCKS_PER_SEC;
   
   if ( run >= 0) runtime_fast_vec.push_back (elapsed_msecs);
   
   //printf (">> computation time: %f ms\n", elapsed_msecs);
   
// ----------------------------------------------------------------------------------------------------

   long double max_abs_error_sqrd = 0.0;
   long double max_rel_error_sqrd = 0.0;
   
   //long double max_abs_weight_error_sqrd = 0.0;
   //long double max_rel_weight_error_sqrd = 0.0;

   int max_abs_error_i = 0;
   
   for ( int i = 0; i < M; ++i )
   {
      //if ( r[i] > 4.0 ) continue;
      
      long double abs_error_sqrd;
      long double rel_error_sqrd;
      
      abs_error_sqrd = ((f_real[i] - f_real_ground_truth[i]) * (f_real[i] - f_real_ground_truth[i]) + (f_imag[i] - f_imag_ground_truth[i]) * (f_imag[i] - f_imag_ground_truth[i]));
      rel_error_sqrd = abs_error_sqrd / (f_real_ground_truth[i] * f_real_ground_truth[i] + f_imag_ground_truth[i] * f_imag_ground_truth[i]);
      
      if ( abs_error_sqrd > max_abs_error_sqrd )
      {
         max_abs_error_sqrd = abs_error_sqrd;
         max_abs_error_i = i;
      }
         
      if ( rel_error_sqrd > max_rel_error_sqrd )
         max_rel_error_sqrd = rel_error_sqrd;
      
      //if ( abs_error_sqrd * exp (- 2.0 * r[i] * r[i]) > max_abs_weight_error_sqrd )
      //   max_abs_weight_error_sqrd = abs_error_sqrd * exp (- 2.0 * r[i] * r[i]);
      //
      //if ( rel_error_sqrd * exp (- 2.0 * r[i] * r[i]) > max_rel_weight_error_sqrd )
      //   max_rel_weight_error_sqrd = rel_error_sqrd * exp (- 2.0 * r[i] * r[i]);      
   }
   
   double max_abs_error = sqrt (max_abs_error_sqrd);
   //double max_abs_weight_error = sqrt (max_abs_weight_error_sqrd);
   double max_rel_error = sqrt (max_rel_error_sqrd);   
   //double max_rel_weight_error = sqrt (max_rel_weight_error_sqrd);
   
   //printf (">> maximum absolute reconstruction error: %.4e \n", max_abs_error);
   //printf ("   (@ r[%d] = %f --- %.15f + i%.15f vs. %.15f + i%.15f)\n", max_abs_error_i, r[max_abs_error_i], f_real_ground_truth[max_abs_error_i], f_imag_ground_truth[max_abs_error_i], f_real[max_abs_error_i], f_imag[max_abs_error_i]);
   
   //printf (">> maximum absolute weighted reconstruction error: %.4e\n", max_abs_weight_error);   
   //printf (">> maximum relative reconstruction error: %.4e\n", max_rel_error);
   //printf (">> maximum relative weighted reconstruction error: %.4e\n\n", max_rel_weight_error);
   
   //printf ("\n");
   
   if ( run >= 0 ) 
   {
      max_abs_error_vec.push_back (max_abs_error);
      max_rel_error_vec.push_back (max_rel_error);
      
      //max_abs_weight_error_vec.push_back (max_abs_weight_error);
      //max_rel_weight_error_vec.push_back (max_rel_weight_error);
   }
   
// ----------------------------------------------------------------------------------------------------   

   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         delete[] f_hat_real[n - 1][l];
         delete[] f_hat_imag[n - 1][l];
      }
      
      delete[] f_hat_real[n - 1];
      delete[] f_hat_imag[n - 1];
   }
   
   delete[] f_hat_real;
   delete[] f_hat_imag;
   
   for( int i = 0; i < M; ++i )
      delete[] x[i];

   delete[] x;
   
   delete[] f_real;
   delete[] f_imag;
   
   delete[] f_real_ground_truth;
   delete[] f_imag_ground_truth;
   
// ----------------------------------------------------------------------------------------------------

   }
   
   double runtime_slow_mean = std::accumulate (runtime_slow_vec.begin (), runtime_slow_vec.end (), 0.0) / runtime_slow_vec.size ();
   double runtime_slow_stdv = sqrt (std::inner_product (runtime_slow_vec.begin (), runtime_slow_vec.end (), runtime_slow_vec.begin (), 0.0) / runtime_slow_vec.size () - runtime_slow_mean * runtime_slow_mean);
   
   printf ("average runtime NDSGLFT: %e +/- %e s\n", runtime_slow_mean / 1000.0, runtime_slow_stdv / 1000.0);   
   
   double runtime_fast_mean = std::accumulate (runtime_fast_vec.begin (), runtime_fast_vec.end (), 0.0) / runtime_fast_vec.size ();
   double runtime_fast_stdv = sqrt (std::inner_product (runtime_fast_vec.begin (), runtime_fast_vec.end (), runtime_fast_vec.begin (), 0.0) / runtime_fast_vec.size () - runtime_fast_mean * runtime_fast_mean);
   
   printf ("average runtime NFSGLFT: %e +/- %e s\n", runtime_fast_mean / 1000.0, runtime_fast_stdv / 1000.0);
   
   double max_abs_error_mean = std::accumulate (max_abs_error_vec.begin (), max_abs_error_vec.end (), 0.0) / max_abs_error_vec.size ();
   double max_abs_error_stdv = sqrt (std::inner_product (max_abs_error_vec.begin (), max_abs_error_vec.end (), max_abs_error_vec.begin (), 0.0) / max_abs_error_vec.size () - max_abs_error_mean * max_abs_error_mean);
   
   printf ("average maximum absolute error NFSGLFT: %e +/- %e\n", max_abs_error_mean, max_abs_error_stdv);
   
   double max_rel_error_mean = std::accumulate (max_rel_error_vec.begin (), max_rel_error_vec.end (), 0.0) / max_rel_error_vec.size ();
   double max_rel_error_stdv = sqrt (std::inner_product (max_rel_error_vec.begin (), max_rel_error_vec.end (), max_rel_error_vec.begin (), 0.0) / max_rel_error_vec.size () - max_rel_error_mean * max_rel_error_mean);
   
   printf ("average maximum relative error NFSGLFT: %e +/- %e\n", max_rel_error_mean, max_rel_error_stdv);
   
   //double max_abs_weight_error_mean = std::accumulate (max_abs_weight_error_vec.begin (), max_abs_weight_error_vec.end (), 0.0) / max_abs_weight_error_vec.size ();
   //double max_abs_weight_error_stdv = sqrt (std::inner_product (max_abs_weight_error_vec.begin (), max_abs_weight_error_vec.end (), max_abs_weight_error_vec.begin (), 0.0) / max_abs_weight_error_vec.size () - max_abs_weight_error_mean * max_abs_weight_error_mean);
   //
   //printf ("average maximum weighted absolute error NFSGLFT-A: %e +/- %e\n", max_abs_weight_error_mean, max_abs_weight_error_stdv);
   //
   //double max_rel_weight_error_mean = std::accumulate (max_rel_weight_error_vec.begin (), max_rel_weight_error_vec.end (), 0.0) / max_rel_weight_error_vec.size ();
   //double max_rel_weight_error_stdv = sqrt (std::inner_product (max_rel_weight_error_vec.begin (), max_rel_weight_error_vec.end (), max_rel_weight_error_vec.begin (), 0.0) / max_rel_weight_error_vec.size () - max_rel_weight_error_mean * max_rel_weight_error_mean);
   //
   //printf ("average maximum weighted relative error NFSGLFT-A: %e +/- %e\n\n", max_rel_weight_error_mean, max_rel_weight_error_stdv);   

   runtime_fast_vec.clear ();
   runtime_slow_vec.clear ();
   
   max_abs_error_vec.clear ();
   max_rel_error_vec.clear ();
   
   //max_abs_weight_error_vec.clear ();
   //max_rel_weight_error_vec.clear ();   
   
   //abs_mean << max_abs_error_mean << ", "; //std::endl;
   //abs_stdv << max_abs_error_stdv << ", "; //std::endl;
   
   //rel_mean << max_rel_error_mean << ", "; //std::endl;
   //rel_stdv << max_rel_error_stdv << ", "; //std::endl;
   
   //abs_mean.flush ();
   //abs_stdv.flush ();
   
   //rel_mean.flush ();
   //rel_stdv.flush ();
   
   }
   
   //abs_mean.close ();
   //abs_stdv.close ();
   
   //rel_mean.close ();
   //rel_stdv.close ();
   
   return 0;
}