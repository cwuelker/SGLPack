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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>

#include "fftw.h"
#include "clenshaw.h"
#include "nfft_use.h"
#include "flt.h"

#include <nfsglft.h>
#include <GSL/gsl-master/specfunc/gsl_sf_legendre.h>

#define NFFT_PRECISION_DOUBLE

using namespace std;

void ndsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag)
{   
   long double R_nl;   
   double* Y_lm = new double[2];
   
   for ( int i = 0; i < M; ++i )
   { 
      f_real[i] = 0.0;
      f_imag[i] = 0.0;
            
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {    
               R_nl = R (n, l, x[i][0]);               
               
               Y (l, m, x[i][1], x[i][2], Y_lm);
                     
               f_real[i] += f_hat_real[n - 1][l][m + l] * R_nl * Y_lm[0];
               f_real[i] -= f_hat_imag[n - 1][l][m + l] * R_nl * Y_lm[1];
                     
               f_imag[i] += f_hat_imag[n - 1][l][m + l] * R_nl * Y_lm[0];
               f_imag[i] += f_hat_real[n - 1][l][m + l] * R_nl * Y_lm[1];
            }
         }
      }
   }

   delete[] Y_lm;
   
   return;
}

void ndsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag)
{   
   long double R_nl;   
   double* Y_lm = new double[2];
   
   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         for ( int m = - l; m <= l; ++m )
         {   
            f_hat_real[n - 1][l][m + l] = 0.0;
            f_hat_imag[n - 1][l][m + l] = 0.0;
   
            for ( int i = 0; i < M; ++i )
            { 
               R_nl = R (n, l, x[i][0]);               
               
               Y (l, m, x[i][1], x[i][2], Y_lm);
                     
               f_hat_real[n - 1][l][m + l] += f_real[i] * R_nl * Y_lm[0];
               f_hat_real[n - 1][l][m + l] += f_imag[i] * R_nl * Y_lm[1];

               f_hat_imag[n - 1][l][m + l] += f_imag[i] * R_nl * Y_lm[0];               
               f_hat_imag[n - 1][l][m + l] -= f_real[i] * R_nl * Y_lm[1];
            }
         }
      }
   }
   
   delete[] Y_lm;
   
   return;
}

void nfsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q)
{      
   const double pi = 3.141592653589793238462643383279502884197169399;
   
   double rho = x[0][0];

   for( int i = 1; i < M; ++i )
      if ( x[i][0] > rho ) rho = x[i][0];

   long double* cos_vec = new long double[2 * B];

   for ( int j = 0; j < 2 * B; ++j )
      cos_vec[j] = 0.5 * rho * (1.0 + cos ((2 * j + 1) * pi / (4 * B)));

   clenshaw function_clenshaw;
   
   long double* data_trans_real = new long double[2 * B];
   long double* data_trans_imag = new long double[2 * B];

   double*** g_real = new double**[B];
   double*** g_imag = new double**[B];

   for ( int l = 0; l < B; ++l )
   {
      g_real[l] = new double*[2 * l + 1];
      g_imag[l] = new double*[2 * l + 1];

      for ( int m = - l; m <= l; ++m )
      {
         g_real[l][m + l] = new double[2 * B];
         g_imag[l][m + l] = new double[2 * B];
      }
   }
   
   for ( int l = 0; l < B; ++l )
   {      
      long double* f_hat_real_clenshaw = new long double[B - l];
      long double* f_hat_imag_clenshaw = new long double[B - l];

      for ( int m = - l; m <= l; ++m )
      {         
         for ( int n = l + 1; n <= B; ++n )
         {
            f_hat_real_clenshaw[n - l - 1] = (long double)f_hat_real[n - 1][l][m + l];
            f_hat_imag_clenshaw[n - l - 1] = (long double)f_hat_imag[n - 1][l][m + l];
         }

         function_clenshaw.DRT_Clenshaw (B, l, cos_vec, f_hat_real_clenshaw, f_hat_imag_clenshaw, data_trans_real, data_trans_imag);

         for ( int j = 0; j < 2 * B; ++j )
         {
            g_real[l][m + l][j] = (double)data_trans_real[j];
            g_imag[l][m + l][j] = (double)data_trans_imag[j];
         }
      }

      delete[] f_hat_real_clenshaw;
      delete[] f_hat_imag_clenshaw;
   }

   delete[] data_trans_real;
   delete[] data_trans_imag;
   
   delete[] cos_vec;

   double*** beta_real = new double**[4 * B];
   double*** beta_imag = new double**[4 * B];
 
   for ( int k = - 2 * B; k < 2 * B; ++k )
   {
      beta_real[k + 2 * B] = new double*[B];
      beta_imag[k + 2 * B] = new double*[B];

      for ( int l = 0; l < B; ++l )
      {
         beta_real[k + 2 * B][l] = new double[2 * l + 1];
         beta_imag[k + 2 * B][l] = new double[2 * l + 1];
      }
   }

   double* f_hat_real_fftw = new double[2 * B];
   double* f_hat_imag_fftw = new double[2 * B];
  
   int index;
   
   fftw function_fftw;
   
   for ( int l = 0; l < B; ++l )
   {      
      for ( int m = - l; m <= l; ++m )
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            f_hat_real_fftw[j] = g_real[l][m + l][j];
            f_hat_imag_fftw[j] = g_imag[l][m + l][j];
         }
              
         function_fftw.dct_mdfd (2 * B, f_hat_real_fftw, f_hat_imag_fftw);
         
         beta_real[0][l][m + l] = 0.0;
         beta_imag[0][l][m + l] = 0.0;

         for ( int k = 1; k < 4 * B; ++k )
         {
            if ( k >= 2 * B )
            {
               index = k - 2 * B;
            }
            else
            {
               index = 2 * B - k;
            }

            if ( k == 2 * B )
            {
               beta_real[k][l][m + l] = f_hat_real_fftw[index];
               beta_imag[k][l][m + l] = f_hat_imag_fftw[index];
            }
            else
            {
               beta_real[k][l][m + l] = 0.5 * f_hat_real_fftw[index];
               beta_imag[k][l][m + l] = 0.5 * f_hat_imag_fftw[index];
            }
         }
      }
   }

   for ( int l = 0; l < B; ++l )
   {
      for ( int m = - l; m <= l; ++m )
      {
         delete[] g_real[l][m + l];
         delete[] g_imag[l][m + l];
      }

      delete[] g_real[l];
      delete[] g_imag[l];
   }

   delete[] g_real;
   delete[] g_imag;   

   double*** h_real = new double**[4 * B];
   double*** h_imag = new double**[4 * B];

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      h_real[k0 + 2 * B] = new double*[2 * B - 1];
      h_imag[k0 + 2 * B] = new double*[2 * B - 1];

      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         h_real[k0 + 2 * B][m + B - 1] = new double[2 * B];
         h_imag[k0 + 2 * B][m + B - 1] = new double[2 * B];
      }
   }

   FLT function_flt;

   double* data_real = new double[2 * B];
   double* data_imag = new double[2 * B];
   
   for ( int m = 1 - B; m <= B - 1; ++m )
   {      
      double* data_FLT_real = new double[B - abs (m)];
      double* data_FLT_imag = new double[B - abs (m)];

      for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
      {         
         for ( int l = abs (m); l < B; ++l )
         {
             data_FLT_real[l - abs (m)] = beta_real[k0 + 2 * B][l][m + l];
             data_FLT_imag[l - abs (m)] = beta_imag[k0 + 2 * B][l][m + l];
         }

         function_flt.flt_seminaive_adjoint (B, m, data_FLT_real, data_real);
         function_flt.flt_seminaive_adjoint (B, m, data_FLT_imag, data_imag);

         for ( int j = 0; j < 2 * B; ++j )
         {
            h_real[k0 + 2 * B][m + B - 1][j] = data_real[j];
            h_imag[k0 + 2 * B][m + B - 1][j] = data_imag[j];
         }
      }

      delete[] data_FLT_real;
      delete[] data_FLT_imag;
   }
   
   delete data_real;
   delete data_imag;

   for ( int k = - 2 * B; k < 2 * B; ++k )
   {
      for ( int l = 0; l < B; ++l )
      {
         delete[] beta_real[k + 2 * B][l];
         delete[] beta_imag[k + 2 * B][l];
      }

      delete[] beta_real[k + 2 * B];
      delete[] beta_imag[k + 2 * B];
   }

   delete[] beta_real;
   delete[] beta_imag;

   double*** zeta_real = new double**[4 * B];
   double*** zeta_imag = new double**[4 * B];

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      zeta_real[k0 + 2 * B] = new double*[4 * B];
      zeta_imag[k0 + 2 * B] = new double*[4 * B];

      for ( int k = - 2 * B; k < 2 * B; ++k )
      {
         zeta_real[k0 + 2 * B][k + 2 * B] = new double[4 * B];
         zeta_imag[k0 + 2 * B][k + 2 * B] = new double[4 * B];

         for (int m = 0; m < 4 * B; ++m )
         {
            zeta_real[k0 + 2 * B][k + 2 * B][m] = 0.0;
            zeta_imag[k0 + 2 * B][k + 2 * B][m] = 0.0;
         }
      }
   }

   double sgnk;

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            f_hat_real_fftw[i] = h_real[k0 + 2 * B][m + B - 1][i];
            f_hat_imag_fftw[i] = h_imag[k0 + 2 * B][m + B - 1][i];
         }
         
         function_fftw.dct_mdfd (2 * B, f_hat_real_fftw, f_hat_imag_fftw);

         zeta_real[k0 + 2 * B][0][m + B - 1] = 0.0;
         zeta_imag[k0 + 2 * B][0][m + B - 1] = 0.0;

         for ( int k = 1; k < 4 * B; ++k )
         {
            if ( k >= 2 * B )
            {
               index = k - 2 * B;
               sgnk = 1.0;
            }
            else
            {
               index = 2 * B - k;
               sgnk = - 1.0;
            }

            if ( m % 2 == 0 )
            {
               if ( k == 2 * B )
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = f_hat_real_fftw[index];
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = f_hat_imag_fftw[index];
               }
               else
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = 0.5 * f_hat_real_fftw[index];
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = 0.5 * f_hat_imag_fftw[index];
               }
            }
            else
            {
               if ( index == 1 )
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = 0.25 * sgnk * (2.0 * f_hat_imag_fftw[0] - f_hat_imag_fftw[2]);
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = - 0.25 * sgnk * (2.0 * f_hat_real_fftw[0] - f_hat_real_fftw[2]);
               }
               else if ( k == 1 || k == 2 || k == 4 * B - 2 || k == 4 * B - 1 )
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = 0.25 * sgnk * f_hat_imag_fftw[index - 1];
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = - 0.25 * sgnk * f_hat_real_fftw[index - 1];
               }
               else if ( k == 2 * B )
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = 0.0;
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = 0.0;
               }
               else
               {
                  zeta_real[k0 + 2 * B][k][m + B - 1] = 0.25 * sgnk * (f_hat_imag_fftw[index - 1] - f_hat_imag_fftw[index + 1]);
                  zeta_imag[k0 + 2 * B][k][m + B - 1] = - 0.25 * sgnk * (f_hat_real_fftw[index - 1] - f_hat_real_fftw[index + 1]);
               }
            }
         }
      }
   }

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         delete[] h_real[k0 + 2 * B][m + B - 1];
         delete[] h_imag[k0 + 2 * B][m + B - 1];
      }

      delete[] h_real[k0 + 2 * B];
      delete[] h_imag[k0 + 2 * B];
   }

   delete[] h_real;
   delete[] h_imag;
   
   delete[] f_hat_real_fftw;
   delete[] f_hat_imag_fftw;
   
   double*** eta_real = new double**[4 * B];
   double*** eta_imag = new double**[4 * B];
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      eta_real[k0 + 2 * B] = new double*[2 * B];
      eta_imag[k0 + 2 * B] = new double*[2 * B];
      
      for ( int k1 = - B; k1 < B; ++k1 )
      {
         eta_real[k0 + 2 * B][k1 + B] = new double[2 * B];
         eta_imag[k0 + 2 * B][k1 + B] = new double[2 * B];
         
         for ( int k2 = - B; k2 < B; ++k2 )
         {
            if ( - B < k2 && k2 < B )
            {
               eta_real[k0 + 2 * B][k1 + B][k2 + B] = zeta_real[k0 + 2 * B][k1 + 2 * B][k2 + B - 1];
               eta_imag[k0 + 2 * B][k1 + B][k2 + B] = zeta_imag[k0 + 2 * B][k1 + 2 * B][k2 + B - 1];
            }
            else
            {
               eta_real[k0 + 2 * B][k1 + B][k2 + B] = 0.0;
               eta_imag[k0 + 2 * B][k1 + B][k2 + B] = 0.0;               
            }
         }
      }
   }
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int m = 0; m < 4 * B; ++m )
      {
         delete[] zeta_real[k0 + 2 * B][m];
         delete[] zeta_imag[k0 + 2 * B][m];
      }

      delete[] zeta_real[k0 + 2 * B];
      delete[] zeta_imag[k0 + 2 * B];
   }

   delete[] zeta_real;
   delete[] zeta_imag;   
   
   complex<double>* nfft_f = new complex<double>[M];

   int N[3];
   
   N[0] = 4 * B;
   N[1] = 2 * B;
   N[2] = 2 * B;

   nfft_use nfft_function;
   
   double** x_trans = new double*[M];
   
   for ( int i = 0; i < M; ++i )
   {
      x_trans[i] = new double[3];
      
      x_trans[i][0] = - acos((2.0 * x[i][0] - rho) / rho) / (2.0 * pi);
      x_trans[i][1] = - x[i][1] / (2.0 * pi);
      x_trans[i][2] = - x[i][2] / (2.0 * pi);

      if ( x_trans[i][0] < - 0.5 ) x_trans[i][0] += 1.0;
      if ( x_trans[i][1] < - 0.5 ) x_trans[i][1] += 1.0;     
      if ( x_trans[i][2] < - 0.5 ) x_trans[i][2] += 1.0;
   }
      
   nfft_function.nfft (N, M, eta_real, eta_imag, x_trans, nfft_f, nfft_sigma, nfft_q);
      
   for ( int i = 0; i < M; ++i ) delete[] x_trans[i];
   
   delete[] x_trans;
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {  
      for ( int k1 = - B; k1 < B; ++k1 )
      {
         delete[] eta_real[k0 + 2 * B][k1 + B];
         delete[] eta_imag[k0 + 2 * B][k1 + B];
      }
      
      delete[] eta_real[k0 + 2 * B];
      delete[] eta_imag[k0 + 2 * B];
   }
   
   delete[] eta_real;
   delete[] eta_imag;
   
   for ( int i = 0; i < M; ++i )
   {
      f_real[i] = nfft_f[i].real ();
      f_imag[i] = nfft_f[i].imag ();
   }
   
   delete[] nfft_f;
   
   return;
}

void nfsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q)
{      
   const double pi = 3.141592653589793238462643383279502884197169399;
   
   double rho = x[0][0];

   for( int i = 1; i < M; ++i )
      if ( x[i][0] > rho ) rho = x[i][0];
   
   int N[3];
   
   N[0] = 4 * B;
   N[1] = 4 * B;
   N[2] = 4 * B;

   nfft_use nfft_function;
   
   double** x_trans = new double*[M];
   
   for ( int i = 0; i < M; ++i )
   {
      x_trans[i] = new double[3];
      
      x_trans[i][0] = - acos((2.0 * x[i][0] - rho) / rho) / (2.0 * pi);
      x_trans[i][1] = - x[i][1] / (2.0 * pi);
      x_trans[i][2] = - x[i][2] / (2.0 * pi);

      if ( x_trans[i][0] < - 0.5 ) x_trans[i][0] += 1.0;
      if ( x_trans[i][1] < - 0.5 ) x_trans[i][1] += 1.0;     
      if ( x_trans[i][2] < - 0.5 ) x_trans[i][2] += 1.0;
   }

   complex<double>* nfft_f = new complex<double>[M];  
   
   for ( int i = 0; i < M; ++i )
   {
      nfft_f[i].real (f_real[i]);
      nfft_f[i].imag (f_imag[i]);      
   }
   
   double*** eta_real = new double**[4 * B];
   double*** eta_imag = new double**[4 * B];

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      eta_real[k0 + 2 * B] = new double*[4 * B];
      eta_imag[k0 + 2 * B] = new double*[4 * B];

      for ( int k = - 2 * B; k < 2 * B; ++k )
      {
         eta_real[k0 + 2 * B][k + 2 * B] = new double[4 * B];
         eta_imag[k0 + 2 * B][k + 2 * B] = new double[4 * B];
      }
   }
   
   nfft_function.nfft_adjoint (N, M, eta_real, eta_imag, x_trans, nfft_f, nfft_sigma, nfft_q);  
   
   for ( int i = 0; i < M; ++i ) delete[] x_trans[i];
   
   delete[] x_trans;
   delete[] nfft_f;
   
   double*** zeta_real = new double**[4 * B];
   double*** zeta_imag = new double**[4 * B];

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      zeta_real[k0 + 2 * B] = new double*[4 * B];
      zeta_imag[k0 + 2 * B] = new double*[4 * B];

      for ( int k = - 2 * B; k < 2 * B; ++k )
      {
         zeta_real[k0 + 2 * B][k + 2 * B] = new double[2 * B - 1];
         zeta_imag[k0 + 2 * B][k + 2 * B] = new double[2 * B - 1];

         for (int m = 1 - B; m <= B - 1; ++m )
         {
            zeta_real[k0 + 2 * B][k + 2 * B][m + B - 1] = eta_real[k0 + 2 * B][k + 2 * B][2 * B + m];
            zeta_imag[k0 + 2 * B][k + 2 * B][m + B - 1] = eta_imag[k0 + 2 * B][k + 2 * B][2 * B + m];
         }
      }
   }   
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int k1 = - 2 * B; k1 < 2 * B; ++k1 )
      {
         delete[] eta_real[k0 + 2 * B][k1 + 2 * B];
         delete[] eta_imag[k0 + 2 * B][k1 + 2 * B];
      }

      delete[] eta_real[k0 + 2 * B];
      delete[] eta_imag[k0 + 2 * B];
   }

   delete[] eta_real;
   delete[] eta_imag;
   
   double*** h_real = new double**[4 * B];
   double*** h_imag = new double**[4 * B];

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      h_real[k0 + 2 * B] = new double*[2 * B - 1];
      h_imag[k0 + 2 * B] = new double*[2 * B - 1];

      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         h_real[k0 + 2 * B][m + B - 1] = new double[2 * B];
         h_imag[k0 + 2 * B][m + B - 1] = new double[2 * B];
      }
   }
   
   double* f_hat_real_fftw = new double[2 * B];
   double* f_hat_imag_fftw = new double[2 * B];
   
   fftw function_fftw;

   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         if ( m % 2 == 0 )
         {
            for ( int i = 0; i < 2 * B; ++i )
            {
               f_hat_real_fftw[i] = 0.5 * (zeta_real[k0 + 2 * B][2 * B - i][m + B - 1] + zeta_real[k0 + 2 * B][2 * B + i][m + B - 1]);
               f_hat_imag_fftw[i] = 0.5 * (zeta_imag[k0 + 2 * B][2 * B - i][m + B - 1] + zeta_imag[k0 + 2 * B][2 * B + i][m + B - 1]);  
            }       
         }
         else
         {            
            f_hat_real_fftw[0] = 0.5 * (zeta_imag[k0 + 2 * B][2 * B - 1][m + B - 1] - zeta_imag[k0 + 2 * B][2 * B + 1][m + B - 1]);
            f_hat_imag_fftw[0] = 0.5 * (zeta_real[k0 + 2 * B][2 * B + 1][m + B - 1] - zeta_real[k0 + 2 * B][2 * B - 1][m + B - 1]);
            
            f_hat_real_fftw[1] = 0.25 * (zeta_imag[k0 + 2 * B][2 * B - 2][m + B - 1] - zeta_imag[k0 + 2 * B][2 * B + 2][m + B - 1]);
            f_hat_imag_fftw[1] = 0.25 * (zeta_real[k0 + 2 * B][2 * B + 2][m + B - 1] - zeta_real[k0 + 2 * B][2 * B - 2][m + B - 1]);            
            
            for ( int i = 2; i < 2 * B - 1; ++i )
            {
               f_hat_real_fftw[i] = 0.25 * (zeta_imag[k0 + 2 * B][2 * B - i - 1][m + B - 1] - zeta_imag[k0 + 2 * B][2 * B - i + 1][m + B - 1]);
               f_hat_real_fftw[i] += 0.25 * (zeta_imag[k0 + 2 * B][2 * B + i - 1][m + B - 1] - zeta_imag[k0 + 2 * B][2 * B + i + 1][m + B - 1]);
               
               f_hat_imag_fftw[i] = 0.25 * (zeta_real[k0 + 2 * B][2 * B - i + 1][m + B - 1] - zeta_real[k0 + 2 * B][2 * B - i - 1][m + B - 1]);
               f_hat_imag_fftw[i] += 0.25 * (zeta_real[k0 + 2 * B][2 * B + i + 1][m + B - 1] - zeta_real[k0 + 2 * B][2 * B + i - 1][m + B - 1]);
            }
            
            f_hat_real_fftw[2 * B - 1] = 0.0;
            f_hat_imag_fftw[2 * B - 1] = 0.0;
         }
         
         function_fftw.dct_mdfd_adjoint (2 * B, f_hat_real_fftw, f_hat_imag_fftw);         
         
         for ( int i = 0; i < 2 * B; ++i )
         {
            h_real[k0 + 2 * B][m + B - 1][i] = f_hat_real_fftw[i];
            h_imag[k0 + 2 * B][m + B - 1][i] = f_hat_imag_fftw[i];
         }         
      }
   }
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int k = - 2 * B; k < 2 * B; ++k )
      {
         delete[] zeta_real[k0 + 2 * B][k + 2 * B];
         delete[] zeta_imag[k0 + 2 * B][k + 2 * B];
      }

      delete[] zeta_real[k0 + 2 * B];
      delete[] zeta_imag[k0 + 2 * B];
   }

   delete[] zeta_real;
   delete[] zeta_imag;

   double*** beta_real = new double**[4 * B];
   double*** beta_imag = new double**[4 * B];
 
   for ( int k = - 2 * B; k < 2 * B; ++k )
   {
      beta_real[k + 2 * B] = new double*[B];
      beta_imag[k + 2 * B] = new double*[B];

      for ( int l = 0; l < B; ++l )
      {
         beta_real[k + 2 * B][l] = new double[2 * l + 1];
         beta_imag[k + 2 * B][l] = new double[2 * l + 1];
      }
   }   
   
   FLT function_flt;

   double* data_real = new double[2 * B];
   double* data_imag = new double[2 * B];
   
   for ( int m = 1 - B; m <= B - 1; ++m )
   {      
      double* data_FLT_real = new double[B - abs (m)];
      double* data_FLT_imag = new double[B - abs (m)];

      for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
      {      
         for ( int j = 0; j < 2 * B; ++j )
         {
            data_real[j] = h_real[k0 + 2 * B][m + B - 1][j];
            data_imag[j] = h_imag[k0 + 2 * B][m + B - 1][j];
         }
            
         function_flt.flt_seminaive (B, m, data_real, data_FLT_real);
         function_flt.flt_seminaive (B, m, data_imag, data_FLT_imag);
         
         for ( int l = abs (m); l < B; ++l )
         {
             beta_real[k0 + 2 * B][l][m + l] = data_FLT_real[l - abs (m)];
             beta_imag[k0 + 2 * B][l][m + l] = data_FLT_imag[l - abs (m)];
         }
      }

      delete[] data_FLT_real;
      delete[] data_FLT_imag;
   }
   
   delete data_real;
   delete data_imag;   
   
   for ( int k0 = - 2 * B; k0 < 2 * B; ++k0 )
   {
      for ( int m = 1 - B; m <= B - 1; ++m )
      {
         delete[] h_real[k0 + 2 * B][m + B - 1];
         delete[] h_imag[k0 + 2 * B][m + B - 1];
      }

      delete[] h_real[k0 + 2 * B];
      delete[] h_imag[k0 + 2 * B];
   }

   delete[] h_real;
   delete[] h_imag;
   
   double*** g_real = new double**[B];
   double*** g_imag = new double**[B];

   for ( int l = 0; l < B; ++l )
   {
      g_real[l] = new double*[2 * l + 1];
      g_imag[l] = new double*[2 * l + 1];

      for ( int m = - l; m <= l; ++m )
      {
         g_real[l][m + l] = new double[2 * B];
         g_imag[l][m + l] = new double[2 * B];
      }
   }   
   
   for ( int l = 0; l < B; ++l )
   {      
      for ( int m = - l; m <= l; ++m )
      {         
         for ( int index = 0; index < 2 * B; ++index )
         {
            f_hat_real_fftw[index] = 0.5 * (beta_real[2 * B - index][l][m + l] + beta_real[2 * B + index][l][m + l]);
            f_hat_imag_fftw[index] = 0.5 * (beta_imag[2 * B - index][l][m + l] + beta_imag[2 * B + index][l][m + l]);
         }
         
         function_fftw.dct_mdfd_adjoint (2 * B, f_hat_real_fftw, f_hat_imag_fftw);         
         
         for ( int j = 0; j < 2 * B; ++j )
         {
            g_real[l][m + l][j] = f_hat_real_fftw[j];
            g_imag[l][m + l][j] = f_hat_imag_fftw[j];
         }
      }
   }   
   
   for ( int k = - 2 * B; k < 2 * B; ++k )
   {
      for ( int l = 0; l < B; ++l )
      {
         delete[] beta_real[k + 2 * B][l];
         delete[] beta_imag[k + 2 * B][l];
      }

      delete[] beta_real[k + 2 * B];
      delete[] beta_imag[k + 2 * B];
   }

   delete[] beta_real;
   delete[] beta_imag;
   
   long double* cos_vec = new long double[2 * B];

   for ( int j = 0; j < 2 * B; ++j )
      cos_vec[j] = 0.5 * rho * (1.0 + cos ((2 * j + 1) * pi / (4 * B)));

   clenshaw function_clenshaw;
   
   long double* data_real_clenshaw = new long double[2 * B];
   long double* data_imag_clenshaw = new long double[2 * B];
   
   for ( int l = 0; l < B; ++l )
   {      
      long double* f_hat_real_clenshaw = new long double[B - l];
      long double* f_hat_imag_clenshaw = new long double[B - l];

      for ( int m = - l; m <= l; ++m )
      {     
         for ( int j = 0; j < 2 * B; ++j )
         {
            data_real_clenshaw[j] = (long double)g_real[l][m + l][j];
            data_imag_clenshaw[j] = (long double)g_imag[l][m + l][j];
         }
         
         function_clenshaw.DRT_Clenshaw_adjoint (B, l, cos_vec, data_real_clenshaw, data_imag_clenshaw, f_hat_real_clenshaw, f_hat_imag_clenshaw);         
         
         for ( int n = l + 1; n <= B; ++n )
         {
            f_hat_real[n - 1][l][m + l] = (double)f_hat_real_clenshaw[n - l - 1];
            f_hat_imag[n - 1][l][m + l] = (double)f_hat_imag_clenshaw[n - l - 1];
         }
      }

      delete[] f_hat_real_clenshaw;
      delete[] f_hat_imag_clenshaw;
   }

   delete[] data_real_clenshaw;
   delete[] data_imag_clenshaw;
   
   delete[] cos_vec;
   
   for ( int l = 0; l < B; ++l )
   {
      for ( int m = - l; m <= l; ++m )
      {
         delete[] g_real[l][m + l];
         delete[] g_imag[l][m + l];
      }

      delete[] g_real[l];
      delete[] g_imag[l];
   }

   delete[] g_real;
   delete[] g_imag;     
   
   delete[] f_hat_real_fftw;
   delete[] f_hat_imag_fftw;   
   
   return;
}

void infsglft_golub_van_loan (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences)
{      
   char* filename_max_abs_error = new char[1000];
   char* filename_max_rel_error = new char[1000];
   char* filename_residue = new char[1000];   
   
   //sprintf (filename_max_abs_error, "/var/run/media/wuelker/CW/Corporate Identity/C - Promotion/Dissertation/graphics/nfsglft_a/trunk/%s_abs", preferences);
   //sprintf (filename_max_rel_error, "/var/run/media/wuelker/CW/Corporate Identity/C - Promotion/Dissertation/graphics/nfsglft_a/trunk/%s_rel", preferences);
   //sprintf (filename_residue, "/var/run/media/wuelker/CW/Corporate Identity/C - Promotion/Dissertation/graphics/nfsglft_a/trunk/%s_res", preferences);   
   
   sprintf (filename_max_abs_error, "/home/wuelker/Desktop/trunk/%s_abs", preferences);
   sprintf (filename_max_rel_error, "/home/wuelker/Desktop/trunk/%s_rel", preferences);
   sprintf (filename_residue, "/home/wuelker/Desktop/trunk/%s_res", preferences);      
   
   std::ofstream max_abs_error_out (filename_max_abs_error, std::ofstream::out);
   std::ofstream max_rel_error_out (filename_max_rel_error, std::ofstream::out);
   std::ofstream residue_out (filename_residue, std::ofstream::out);   
   
   max_abs_error_out << std::scientific << std::setprecision (10);
   max_rel_error_out << std::scientific << std::setprecision (10);
   residue_out << std::scientific << std::setprecision (10);   
   
   delete[] filename_max_abs_error;
   delete[] filename_max_rel_error;
   delete[] filename_residue;   
   
   double* r_real = new double[M];
   double* r_imag = new double[M];
   
   double* r_real_temp = new double[M];
   double* r_imag_temp = new double[M];
   
   double* A_p_real = new double[M];
   double* A_p_imag = new double[M];   
   
   double*** p_real = new double**[B];
   double*** p_imag = new double**[B];   
   
   double*** A_T_r_real = new double**[B];
   double*** A_T_r_imag = new double**[B];   
   
   for ( int n = 1; n <= B; ++n )
   {
      p_real[n - 1] = new double*[n];
      p_imag[n - 1] = new double*[n];
      
      A_T_r_real[n - 1] = new double*[n];
      A_T_r_imag[n - 1] = new double*[n];       
   
      for ( int l = 0; l < n; ++l )
      {
         p_real[n - 1][l] = new double[2 * l + 1];
         p_imag[n - 1][l] = new double[2 * l + 1];
         
         A_T_r_real[n - 1][l] = new double[2 * l + 1];
         A_T_r_imag[n - 1][l] = new double[2 * l + 1];
      }
   }
   
   nfsglft (B, M, x, r_real, r_imag, f_hat_real, f_hat_imag, nfft_sigma, nfft_q);
   
   long double alpha, beta, enumerator, denominator;
   
   for ( int i = 0; i < M; ++i )
   {
      r_real[i] = f_real[i] - r_real[i];
      r_imag[i] = f_imag[i] - r_imag[i];
   }
   
   for ( int iter = 0; iter < n_iter; ++iter )
   {
      printf ("\n%d/%d... ", iter + 1, n_iter);
      fflush(stdout);
      
      if ( iter == 0 )
      {
         nfsglft_adjoint (B, M, x, r_real, r_imag, p_real, p_imag, nfft_sigma, nfft_q);
         
         enumerator = 0.0;
         
         for ( int n = 1; n <= B; ++n )
         {
            for ( int l = 0; l < n; ++l )
            {
               for ( int m = - l; m <= l; ++m )
               {
                  enumerator += (long double)(p_real[n - 1][l][m + l] * p_real[n - 1][l][m + l] + p_imag[n - 1][l][m + l] * p_imag[n - 1][l][m + l]);
               }
            }
         }
      }
      else
      {                  
         denominator = enumerator;
         
         nfsglft_adjoint (B, M, x, r_real, r_imag, A_T_r_real, A_T_r_imag, nfft_sigma, nfft_q);
         
         enumerator = 0.0;
         
         for ( int n = 1; n <= B; ++n )
         {
            for ( int l = 0; l < n; ++l )
            {
               for ( int m = - l; m <= l; ++m )
               {
                  enumerator += (long double)(A_T_r_real[n - 1][l][m + l] * A_T_r_real[n - 1][l][m + l] + A_T_r_imag[n - 1][l][m + l] * A_T_r_imag[n - 1][l][m + l]);
               }
            }
         }

         beta = enumerator / denominator;
         
         for ( int n = 1; n <= B; ++n )
         {
            for ( int l = 0; l < n; ++l )
            {
               for ( int m = - l; m <= l; ++m )
               {
                  p_real[n - 1][l][m + l] *= beta;
                  p_imag[n - 1][l][m + l] *= beta;
                  
                  p_real[n - 1][l][m + l] += A_T_r_real[n - 1][l][m + l];
                  p_imag[n - 1][l][m + l] += A_T_r_imag[n - 1][l][m + l];                  
               }
            }
         }
      }
      
      nfsglft (B, M, x, A_p_real, A_p_imag, p_real, p_imag, nfft_sigma, nfft_q);
      
      denominator = 0.0;
      
      for ( int i = 0; i < M; ++i )
         denominator += (long double)(A_p_real[i] * A_p_real[i] + A_p_imag[i] * A_p_imag[i]);
      
      alpha = enumerator / denominator;
      
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {
               f_hat_real[n - 1][l][m + l] += alpha * p_real[n - 1][l][m + l];
               f_hat_imag[n - 1][l][m + l] += alpha * p_imag[n - 1][l][m + l];
            }
         }
      }
      
      for ( int i = 0; i < M; ++i )
      {
         r_real_temp[i] = r_real[i];
         r_imag_temp[i] = r_imag[i];
         
         r_real[i] -= alpha * A_p_real[i];
         r_imag[i] -= alpha * A_p_imag[i];
      }
      
      long double max_abs_error_sqrd = 0.0;
      long double max_rel_error_sqrd = 0.0;
   
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {            
               long double abs_error_sqrd, rel_error_sqrd;

               abs_error_sqrd = ((f_hat_real[n - 1][l][m + l] - f_hat_real_ground_truth[n - 1][l][m + l]) * (f_hat_real[n - 1][l][m + l] - f_hat_real_ground_truth[n - 1][l][m + l]) + (f_hat_imag[n - 1][l][m + l] - f_hat_imag_ground_truth[n - 1][l][m + l]) * (f_hat_imag[n - 1][l][m + l] - f_hat_imag_ground_truth[n - 1][l][m + l]));
               rel_error_sqrd = abs_error_sqrd / (f_hat_real_ground_truth[n - 1][l][m + l] * f_hat_real_ground_truth[n - 1][l][m + l] + f_hat_imag_ground_truth[n - 1][l][m + l] * f_hat_imag_ground_truth[n - 1][l][m + l]);            
            
               if ( abs_error_sqrd > max_abs_error_sqrd ) max_abs_error_sqrd = abs_error_sqrd;
               if ( rel_error_sqrd > max_rel_error_sqrd ) max_rel_error_sqrd = rel_error_sqrd;
            }
         }
      }
   
      double max_abs_error = sqrt (max_abs_error_sqrd);
      double max_rel_error = sqrt (max_rel_error_sqrd);   

      max_abs_error_out << max_abs_error << std::endl;
      max_rel_error_out << max_rel_error << std::endl;
      
      max_abs_error_out.flush ();
      max_rel_error_out.flush ();
      
      double norm_r_sqrd = 0.0;
      
      for ( int i = 0; i < M; ++i )
         norm_r_sqrd += r_real[i] * r_real[i] + r_imag[i] * r_imag[i];
      
      double residue = sqrt (norm_r_sqrd);
      
      residue_out << residue << std::endl;
      
      residue_out.flush ();
      
      printf ("(residuum: %.2e, max. abs. / rel. error: %.2e / %.2e, alpha / beta: %.2Le / %.2Le)", residue, max_abs_error, max_rel_error, alpha, beta);
   }
   
   delete[] r_real;
   delete[] r_imag;
   
   delete[] r_real_temp;
   delete[] r_imag_temp;
   
   delete[] A_p_real;
   delete[] A_p_imag;   
   
   for ( int n = 1; n <= B; ++n )
   {
      for ( int l = 0; l < n; ++l )
      {
         delete[] p_real[n - 1][l];
         delete[] p_imag[n - 1][l];         
         
         delete[] A_T_r_real[n - 1][l];
         delete[] A_T_r_imag[n - 1][l];  
      }
      
      delete[] p_real[n - 1];
      delete[] p_imag[n - 1];      
      
      delete[] A_T_r_real[n - 1];
      delete[] A_T_r_imag[n - 1]; 
   }
   
   delete[] p_real;
   delete[] p_imag;   
   
   delete[] A_T_r_real;
   delete[] A_T_r_imag;
   
   max_abs_error_out.close ();
   max_rel_error_out.close ();
   residue_out.close ();   
   
   return;
}

void infsglft_bjoerck (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences)
{      
   double* r_k_real = new double[M];
   double* r_k_imag = new double[M];
   
   ndsglft (B, M, x, r_k_real, r_k_imag, f_hat_real, f_hat_imag);
   
   for ( int i = 0; i < M; ++i )
   {
      r_k_real[i] = f_real[i] - r_k_real[i];
      r_k_imag[i] = f_imag[i] - r_k_imag[i];
   }
   
   double*** p_k_real = new double**[B];
   double*** p_k_imag = new double**[B];     
   
   for ( int n = 1; n <= B; ++n )
   {
      p_k_real[n - 1] = new double*[n];
      p_k_imag[n - 1] = new double*[n];   
   
      for ( int l = 0; l < n; ++l )
      {
         p_k_real[n - 1][l] = new double[2 * l + 1];
         p_k_imag[n - 1][l] = new double[2 * l + 1];       
      }
   }
   
   ndsglft_adjoint (B, M, x, r_k_real, r_k_imag, p_k_real, p_k_imag);

   double*** s_k_real = new double**[B];
   double*** s_k_imag = new double**[B];

   double gamma_k = 0.0;
   
   for ( int n = 1; n <= B; ++n )
   {
      s_k_real[n - 1] = new double*[n];
      s_k_imag[n - 1] = new double*[n];   
   
      for ( int l = 0; l < n; ++l )
      {
         s_k_real[n - 1][l] = new double[2 * l + 1];
         s_k_imag[n - 1][l] = new double[2 * l + 1];
         
         for ( int m = - l; m <= l; ++m )
         {
            s_k_real[n - 1][l][m + l] = p_k_real[n - 1][l][m + l];
            s_k_imag[n - 1][l][m + l] = p_k_imag[n - 1][l][m + l];
            
            gamma_k += s_k_real[n - 1][l][m + l] * s_k_real[n - 1][l][m + l] + s_k_imag[n - 1][l][m + l] * s_k_imag[n - 1][l][m + l];
         }
      }
   }
   
   double beta_k, gamma_k_plus_1;
   
   double* q_k_real = new double[M];
   double* q_k_imag = new double[M];
   
   for ( int k = 0; k < n_iter; ++k )
   {      
      printf ("\n%d/%d... ", k + 1, n_iter);
      fflush(stdout);
      
      ndsglft (B, M, x, q_k_real, q_k_imag, p_k_real, p_k_imag);
      
      double alpha_k = 0.0;
      
      for ( int i = 0; i < M; ++i ) alpha_k += q_k_real[i] * q_k_real[i] + q_k_imag[i] * q_k_imag[i];
      
      alpha_k = gamma_k / alpha_k;
      
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {  
            for ( int m = - l; m <= l; ++m )
            {
               f_hat_real[n - 1][l][m + l] += alpha_k * p_k_real[n - 1][l][m + l];
               f_hat_imag[n - 1][l][m + l] += alpha_k * p_k_imag[n - 1][l][m + l];
            }
         }
      }
      
      for ( int i = 0; i < M; ++i )
      {
         r_k_real[i] -= alpha_k * q_k_real[i];
         r_k_imag[i] -= alpha_k * q_k_imag[i];
      }
      
      ndsglft_adjoint (B, M, x, r_k_real, r_k_imag, s_k_real, s_k_imag);
      
      gamma_k_plus_1 = 0.0;
      
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {       
               gamma_k_plus_1 += s_k_real[n - 1][l][m + l] * s_k_real[n - 1][l][m + l] + s_k_imag[n - 1][l][m + l] * s_k_imag[n - 1][l][m + l];
            }
         }
      }
      
      beta_k = gamma_k_plus_1 / gamma_k;
      
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {       
               p_k_real[n - 1][l][m + l] = s_k_real[n - 1][l][m + l] + beta_k * p_k_real[n - 1][l][m + l];
               p_k_imag[n - 1][l][m + l] = s_k_imag[n - 1][l][m + l] + beta_k * p_k_imag[n - 1][l][m + l];
            }
         }
      }
      
      gamma_k = gamma_k_plus_1;
      
// ----------------------------------------------------------------------------------------------------

      long double max_abs_error_sqrd = 0.0;
      long double max_rel_error_sqrd = 0.0;
   
      for ( int n = 1; n <= B; ++n )
      {
         for ( int l = 0; l < n; ++l )
         {
            for ( int m = - l; m <= l; ++m )
            {            
               long double abs_error_sqrd, rel_error_sqrd;

               abs_error_sqrd = ((f_hat_real[n - 1][l][m + l] - f_hat_real_ground_truth[n - 1][l][m + l]) * (f_hat_real[n - 1][l][m + l] - f_hat_real_ground_truth[n - 1][l][m + l]) + (f_hat_imag[n - 1][l][m + l] - f_hat_imag_ground_truth[n - 1][l][m + l]) * (f_hat_imag[n - 1][l][m + l] - f_hat_imag_ground_truth[n - 1][l][m + l]));
               rel_error_sqrd = abs_error_sqrd / (f_hat_real_ground_truth[n - 1][l][m + l] * f_hat_real_ground_truth[n - 1][l][m + l] + f_hat_imag_ground_truth[n - 1][l][m + l] * f_hat_imag_ground_truth[n - 1][l][m + l]);            
            
               if ( abs_error_sqrd > max_abs_error_sqrd ) max_abs_error_sqrd = abs_error_sqrd;
               if ( rel_error_sqrd > max_rel_error_sqrd ) max_rel_error_sqrd = rel_error_sqrd;
            }
         }
      }
   
      double max_abs_error = sqrt (max_abs_error_sqrd);
      double max_rel_error = sqrt (max_rel_error_sqrd);
      
      double norm_r_k_sqrd = 0.0;
      
      for ( int i = 0; i < M; ++i ) norm_r_k_sqrd += r_k_real[i] * r_k_real[i] + r_k_imag[i] * r_k_imag[i];
      
      double residue = sqrt (norm_r_k_sqrd);
      
      printf ("(residuum: %.2e, max. abs. / rel. error: %.2e / %.2e)", residue, max_abs_error, max_rel_error);

// ----------------------------------------------------------------------------------------------------      
   }
   
   delete[] q_k_real;
   delete[] q_k_imag;
   
   for ( int n = 1; n <= B; ++n )
   {   
      for ( int l = 0; l < n; ++l )
      {
         delete[] p_k_real[n - 1][l];
         delete[] p_k_imag[n - 1][l];     
         
         delete[] s_k_real[n - 1][l];
         delete[] s_k_imag[n - 1][l]; 
      }
      
      delete[] p_k_real[n - 1];
      delete[] p_k_imag[n - 1]; 
      
      delete[] s_k_real[n - 1];
      delete[] s_k_imag[n - 1]; 
   }
   
   delete[] p_k_real;
   delete[] p_k_imag;
   
   delete[] s_k_real;
   delete[] s_k_imag;   
   
   return;
}

long double R (int n, int l, long double r)
{
   const long double sqrt_sqrt_pi = 1.33133536380038971279753491795;
   
   int n_l_difference = n - l;
   
   long double r_squared = r * r;
   
   long double R;
   
   long double temp[n_l_difference];
   
   temp[0] = 1.0;
   
   if ( n_l_difference > 1 )
   {
      temp[1] = 1.5 + l - r_squared;
      
      for ( int index = 2; index < n_l_difference; ++index )
      {
         temp[index] = 0.0;
      }
      
      for ( int index = 2; index < n_l_difference; ++index )
      {
         temp[index] = (((2.0 * index - 0.5 + l - r_squared) * temp[index - 1]) + ((0.5 - index - l) * temp[index - 2])) / index;
      }
   }
   
   R = temp[n_l_difference - 1];
   
   //R *= pow (r, l);
   for ( int index = 0; index < l; ++index )
      R *= r;
   
   long double factor = pow (2.0, 0.5 * n + 0.5) / sqrt_sqrt_pi;
   
   for ( int index = 2 * n - 1; index > 1; index -= 2 )
   {
      factor /= sqrt (index); 
   }
   
   for ( int index = 2; index < n_l_difference; ++index )
   {
      factor *= sqrt (index); 
   }
   
   R *= factor;

   return R;
}

double M_P (int l, int m, double xi)
{
   double M_lm_P_lm_xi;
   
   if ( m >= 0 )
   {
      M_lm_P_lm_xi = gsl_sf_legendre_sphPlm (l, m, xi);
   }
   else
   {
      M_lm_P_lm_xi = gsl_sf_legendre_sphPlm (l, - m, xi);
      
      if ( - m % 2 )
         M_lm_P_lm_xi = - M_lm_P_lm_xi;
   }
        
   return M_lm_P_lm_xi;
}

void Y (int l, int m, double theta, double phi, double* Y_lm)
{        
   double M_lm_P_lm_cos_theta = M_P (l, m, cos (theta));     
        
   Y_lm[0] = M_lm_P_lm_cos_theta * cos (m * phi);
   Y_lm[1] = M_lm_P_lm_cos_theta * sin (m * phi);
}
