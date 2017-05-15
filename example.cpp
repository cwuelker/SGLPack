//
//        SGLPack 1.0 - Fast spherical Gauss-Laguerre Fourier transforms        
//   
//  
//   Contact: Christian Wuelker, M.Sc.
//            wuelker@math.uni-luebeck.de
//  
//   Copyright 2017  Christian Wuelker, Juergen Prestin
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

#include <fsglft.h>
#include <fsft.h>

int main ()
{
   const int B = 32;

   std::string path = ".";
   
   printf ("\n//----------------------------------------------------------------//");
   printf ("\n// SGLPack 1.0 - Fast spherical Gauss-Laguerre Fourier transforms //");
   printf ("\n// Copyright 2017 Christian Wuelker, Juergen Prestin              //");
   printf ("\n//----------------------------------------------------------------//");
   
   printf ("\n\n>> Bandwidth B = %d\n", B);
        
   double*** f_real = new double**[2 * B];
   double*** f_imag = new double**[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      f_real[i] = new double*[2 * B];
      f_imag[i] = new double*[2 * B];
      
      for ( int j = 0; j < 2 * B; ++j )
      {
         f_real[i][j] = new double[2 * B];
         f_imag[i][j] = new double[2 * B];
      }
   }
   
   double*** f_hat_real = new double**[B];
   double*** f_hat_imag = new double**[B];
   
   double*** f_hat_recon_real = new double**[B];
   double*** f_hat_recon_imag = new double**[B];
   
   for ( int n = 1; n <= B; ++n )
   {
      f_hat_real[n - 1] = new double*[n];
      f_hat_imag[n - 1] = new double*[n];
      
      f_hat_recon_real[n - 1] = new double*[n];
      f_hat_recon_imag[n - 1] = new double*[n];
   
      for ( int l = 0; l < n; ++l )
      {
         f_hat_real[n - 1][l] = new double[2 * l + 1];
         f_hat_imag[n - 1][l] = new double[2 * l + 1];
         
         f_hat_recon_real[n - 1][l] = new double[2 * l + 1];
         f_hat_recon_imag[n - 1][l] = new double[2 * l + 1];
      }
   }
   
// Generate (random) SGL Fourier coefficients
// ----------------------------------------------------------------------------------------------------

   printf (">> Generating random SGL Fourier coefficients... ");

   srand (time (NULL));

   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         for ( int m = - l; m <= l; ++m )
         {
            f_hat_real[n - 1][l][m + l] = (double)(((rand () % 200000) - 100000) / 100000.0);
            f_hat_imag[n - 1][l][m + l] = (double)(((rand () % 200000) - 100000) / 100000.0);
         }
      }
   }
   
   printf ("done.\n");

// ----------------------------------------------------------------------------------------------------

   clock_t begin = clock ();
   
// ----------------------------------------------------------------------------------------------------      
   
   printf (">> Reconstructing function values from SGL Fourier coefficients using the iFSGLFT... ");
   
   ifsglft (B, f_real, f_imag, f_hat_real, f_hat_imag, path);
   
   printf ("done.\n");
   
   printf (">> Reconstructing SGL Fourier coeffcients from reconstructed function values using the FSGLFT... ");
   
   fsglft (B, f_real, f_imag, f_hat_recon_real, f_hat_recon_imag, path);
   
   printf ("done.\n");
   
// ----------------------------------------------------------------------------------------------------   
   
   clock_t end = clock ();
   
// ----------------------------------------------------------------------------------------------------   
   
   double error_sqrd, abs_f_hat_sqrd;
   double abs_error_sqrd = 0.0;
   double rel_error_sqrd = 0.0;
   
   int n_temp, n_temp_2;
   int l_temp, l_temp_2;
   int m_temp, m_temp_2;
   
   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         for ( int m = - l; m <= l; ++m )
         {            
            error_sqrd = pow (f_hat_real[n - 1][l][m + l] - f_hat_recon_real[n - 1][l][m + l], 2.0) + pow (f_hat_imag[n - 1][l][m + l] - f_hat_recon_imag[n - 1][l][m + l], 2.0);
            
            if ( error_sqrd > abs_error_sqrd )
            {
               abs_error_sqrd = error_sqrd;
                
               n_temp = n;
               l_temp = l;
               m_temp = m;
            }
             
            abs_f_hat_sqrd = pow (f_hat_real[n - 1][l][m + l], 2.0) + pow (f_hat_imag[n - 1][l][m + l], 2.0);
             
            if ( error_sqrd / abs_f_hat_sqrd > rel_error_sqrd )
            {
               rel_error_sqrd = error_sqrd / abs_f_hat_sqrd;
                
               n_temp_2 = n;
               l_temp_2 = l;
               m_temp_2 = m;
            }
         }
      }
   }
   
   double abs_error = sqrt (abs_error_sqrd);
   double rel_error = sqrt (rel_error_sqrd);
   
   printf (">> Maximum absolute reconstruction error: %e (at n = %d, l = %d, m = %d)\n", abs_error, n_temp, l_temp, m_temp);
   printf (">> Maximum relative reconstruction error: %e (at n = %d, l = %d, m = %d)\n", rel_error, n_temp_2, l_temp_2, m_temp_2);
   
   for ( int i = 0; i < 2 * B; ++i )
   {      
      for ( int j = 0; j < 2 * B; ++j )
      {
         delete[] f_real[i][j];
         delete[] f_imag[i][j];
      }
      
      delete[] f_real[i];
      delete[] f_imag[i];
   }
   
   delete[] f_real;
   delete[] f_imag;
   
   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         delete[] f_hat_real[n - 1][l];
         delete[] f_hat_imag[n - 1][l];
         
         delete[] f_hat_recon_real[n - 1][l];
         delete[] f_hat_recon_imag[n - 1][l];
      }
      
      delete[] f_hat_real[n - 1];
      delete[] f_hat_imag[n - 1];
      
      delete[] f_hat_recon_real[n - 1];
      delete[] f_hat_recon_imag[n - 1];
   }
   
   delete[] f_hat_real;
   delete[] f_hat_imag;
      
   delete[] f_hat_recon_real;
   delete[] f_hat_recon_imag;   
   
   double elapsed_msecs = 1000.0 * (double)(end - begin) / CLOCKS_PER_SEC;
   
   printf (">> Computation time: %f ms\n\n", elapsed_msecs);
   
   return 0;
}
