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

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string.h>
#include <sstream>

#include <fsft.h>
#include <fsglft.h>

void fsglft (int B, double*** f_real, double*** f_imag, double*** f_hat_real, double*** f_hat_imag, std::string path)
{
   double*** f_fsft_real = new double**[2 * B];
   double*** f_fsft_imag = new double**[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      f_fsft_real[i] = new double*[B];
      f_fsft_imag[i] = new double*[B];
   
      for ( int l = 0; l < B; ++l )
      {
         f_fsft_real[i][l] = new double[2 * B + 1]; 
         f_fsft_imag[i][l] = new double[2 * B + 1];
      }
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      fsft (B, f_real[i], f_imag[i], f_fsft_real[i], f_fsft_imag[i]);
   }
   
   char* filename_r = new char[1000];   
   char* filename_a = new char[1000];
   
   sprintf (filename_r, "%s/precomp/half_range_Hermite_quadrature/%d_100_digits_abcissae.txt", path.c_str (), 2 * B);
   sprintf (filename_a, "%s/precomp/half_range_Hermite_quadrature/%d_1000_digits_weights.txt", path.c_str (), 2 * B);
   
   long double* r = new long double[2 * B];
   long double* a = new long double[2 * B];
   
   std::ifstream in_r (filename_r);
   std::ifstream in_a (filename_a);
   
   std::string line_r;
   std::string line_a;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      std::getline (in_r, line_r);
      std::getline (in_a, line_a);
      
      std::stringstream str_stream_r (line_r);
      std::stringstream str_stream_a (line_a);
         
      str_stream_r >> r[i];
      str_stream_a >> a[i];
   }
   
   delete[] filename_r;
   delete[] filename_a;
   
   long double* rdata_Clenshaw = new long double[2 * B];
   long double* idata_Clenshaw = new long double[2 * B];
   
   long double* data_trans_real;
   long double* data_trans_imag;
   
   for ( int l = 0; l < B; ++l )
   {
      data_trans_real = new long double[B - l];
      data_trans_imag = new long double[B - l];
   
      for ( int m = - l; m <= l; ++m )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            rdata_Clenshaw[i] = a[i] * r[i] * r[i] * f_fsft_real[i][l][m + B - 1];
            idata_Clenshaw[i] = a[i] * r[i] * r[i] * f_fsft_imag[i][l][m + B - 1];
         }
      
         DRT_Clenshaw (B, l, r, rdata_Clenshaw, idata_Clenshaw, data_trans_real, data_trans_imag); 
      
         for ( int n = l + 1; n <= B; ++n )
         {
            f_hat_real[n - 1][l][m + l] = data_trans_real[n - l - 1];
            f_hat_imag[n - 1][l][m + l] = data_trans_imag[n - l - 1];
         }
      }
   
      delete[] data_trans_real;
      delete[] data_trans_imag;
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {      
      for ( int l = 0; l < B; ++l )
      {
         delete[] f_fsft_real[i][l]; 
         delete[] f_fsft_imag[i][l];
      }
      
      delete[] f_fsft_real[i];
      delete[] f_fsft_imag[i];
   }
   
   delete[] f_fsft_real;
   delete[] f_fsft_imag;
   
   delete[] rdata_Clenshaw;
   delete[] idata_Clenshaw;
   delete[] r;
   delete[] a;
   
   return;
}

void ifsglft (int B, double*** f_real, double*** f_imag, double*** f_hat_real, double*** f_hat_imag, std::string path)
{
   double*** f_fsft_real = new double**[2 * B];
   double*** f_fsft_imag = new double**[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      f_fsft_real[i] = new double*[2 * B];
      f_fsft_imag[i] = new double*[2 * B];
   
      for ( int j = 0; j < 2 * B; ++j )
      {
         f_fsft_real[i][j] = new double[2 * B]; 
         f_fsft_imag[i][j] = new double[2 * B];
      }
   }
   
   char* filename_r = new char[1000];   
   
   sprintf (filename_r, "%s/precomp/half_range_Hermite_quadrature/%d_100_digits_abcissae.txt", path.c_str (), 2 * B);
   
   long double* r = new long double[2 * B];
   
   std::ifstream in_r (filename_r);
   
   std::string line_r;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      std::getline (in_r, line_r);
      
      std::stringstream str_stream_r (line_r);
         
      str_stream_r >> r[i];
   }
   
   delete[] filename_r;
   
   long double* data_trans_real = new long double[2 * B];
   long double* data_trans_imag = new long double[2 * B];
   
   long double* rdata_Clenshaw;
   long double* idata_Clenshaw;
   
   for ( int l = 0; l < B; ++l )
   {
      rdata_Clenshaw = new long double[B - l];
      idata_Clenshaw = new long double[B - l];
      
      for ( int m = - l; m <= l; ++m )
      {
         for ( int n = l + 1; n <= B; ++n )
         {
            rdata_Clenshaw[n - l - 1] = f_hat_real[n - 1][l][m + l];
            idata_Clenshaw[n - l - 1] = f_hat_imag[n - 1][l][m + l];
         }
         
         iDRT_Clenshaw (B, l, r, rdata_Clenshaw, idata_Clenshaw, data_trans_real, data_trans_imag);
         
         for ( int i = 0; i < 2 * B; ++i )
         {
            f_fsft_real[i][l][m + l] = data_trans_real[i];
            f_fsft_imag[i][l][m + l] = data_trans_imag[i];
         }
      }
      
      delete[] rdata_Clenshaw;
      delete[] idata_Clenshaw;
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {   
      ifsft (B, f_real[i], f_imag[i], f_fsft_real[i], f_fsft_imag[i]);
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {      
      for ( int j = 0; j < 2 * B; ++j )
      {
         delete[] f_fsft_real[i][j]; 
         delete[] f_fsft_imag[i][j];
      }
      
      delete[] f_fsft_real[i];
      delete[] f_fsft_imag[i];
   }
   
   delete[] f_fsft_real;
   delete[] f_fsft_imag;   
   
   delete[] r;
   delete[] data_trans_real;
   delete[] data_trans_imag;
   
   return;
}

long double R (int n, int l, long double r)
{
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
   
   long double factor = pow (2.0, 0.5 * n + 0.5) / (long double)sqrt(1.77245385090551588191);
   
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

long double R_l_plus_1_l (int l, long double r)
{
   return R (l + 1, l, r);
}

long double R_l_plus_2_l (int l, long double r)
{
   return R (l + 2, l, r);
}

long double alpha (int n, int l, long double r_squared)
{
   return (2 * n - l - 0.5 - r_squared) / (long double)sqrt ((n + 0.5) * (n - l));
}

long double beta (int n, int l)
{
   return - (long double)sqrt ((long double)(n - 0.5) * (long double)(n - l - 1) / (long double)((n + 0.5) * (n - l)));
}

void DRT_Clenshaw (int B, int l, long double* r, long double* rdata, long double* idata, long double* data_trans_real, long double* data_trans_imag)
{
   long double* r_squared = new long double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      r_squared[i] = r[i] * r[i];
      
   long double* data_temp_real = new long double[4 * B];
   long double* data_temp_imag = new long double[4 * B];

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[2 * B + i] = R_l_plus_1_l (l, r[i]);
      data_temp_imag[2 * B + i] = data_temp_real[2 * B + i];
   }
      
   data_trans_real[0] = 0.0;
   data_trans_imag[0] = 0.0;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_trans_real[0] += data_temp_real[2 * B + i] * rdata[i];
      data_trans_imag[0] += data_temp_real[2 * B + i] * idata[i];
   }
      
   if ( B - l == 1 )
      return;

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[i] = R_l_plus_2_l (l, r[i]);
      data_temp_imag[i] = data_temp_real[i];
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[i] = rdata[i] * data_temp_real[i];
      data_temp_imag[i] = idata[i] * data_temp_imag[i];
   }
  
   data_trans_real[1] = 0.0;
   data_trans_imag[1] = 0.0;   
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_trans_real[1] += data_temp_real[i];
      data_trans_imag[1] += data_temp_imag[i];
   }
      
   if ( B - l == 2 )
      return;

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[2 * B + i] = beta (l + 2, l) * rdata[i] * data_temp_real[2 * B + i];
      data_temp_imag[2 * B + i] = beta (l + 2, l) * idata[i] * data_temp_imag[2 * B + i];
   }
   
   long double* temp_real = new long double[2 * B];
   long double* temp_imag = new long double[2 * B];
   
   for ( int index = 1; index <= B - l - 2; ++index )
   {
      if ( index < B - l - 2 )   
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            temp_real[i] = data_temp_real[i];
            temp_imag[i] = data_temp_imag[i];
         }
      }

      for ( int i = 0; i < 2 * B; ++i )
      {
         data_temp_real[i] = alpha (l + 1 + index, l, r_squared[i]) * data_temp_real[i] + data_temp_real[2 * B + i];
         data_temp_imag[i] = alpha (l + 1 + index, l, r_squared[i]) * data_temp_imag[i] + data_temp_imag[2 * B + i];
      }
   
      data_trans_real[2 + index - 1] = 0.0;
      data_trans_imag[2 + index - 1] = 0.0;
      
      for ( int i = 0; i < 2 * B; ++i )
      {
         data_trans_real[2 + index - 1] += data_temp_real[i];
         data_trans_imag[2 + index - 1] += data_temp_imag[i];
      }
      
      if ( index < B - l - 2 )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            data_temp_real[2 * B + i] = beta (l + 2 + index, l) * temp_real[i];
            data_temp_imag[2 * B + i] = beta (l + 2 + index, l) * temp_imag[i];
         }
      }
   }
   
   delete[] data_temp_real;
   delete[] data_temp_imag;
   delete[] temp_real;
   delete[] temp_imag;
   delete[] r_squared;
   
   return;
}

void iDRT_Clenshaw (int B, int l, long double* r, long double* rdata, long double* idata, long double* data_trans_real, long double* data_trans_imag)
{
   long double* r_squared = new long double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      r_squared[i] = r[i] * r[i];
   
   long double* b_plus_1_real = new long double[2 * B];
   long double* b_plus_1_imag = new long double[2 * B];
   
   long double* b_plus_2_real = new long double[2 * B];
   long double* b_plus_2_imag = new long double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      b_plus_1_real[i] = 0.0;
      b_plus_1_imag[i] = 0.0;
      
      b_plus_2_real[i] = 0.0;
      b_plus_2_imag[i] = 0.0;
   }
   
   long double* b_temp_real;
   long double* b_temp_imag;
   
   if ( B - l - 2 >= 1 )
   {
      b_temp_real = new long double[2 * B];
      b_temp_imag = new long double[2 * B];
   }
      
   for ( int index = 0; index < B - l - 1; ++index )
   {
      if ( index > 0 )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            b_temp_real[i] = b_plus_1_real[i];
            b_temp_imag[i] = b_plus_1_imag[i];
            
            b_plus_1_real[i] = alpha (B - index, l, r_squared[i]) * b_plus_1_real[i] + b_plus_2_real[i] + rdata[B - l - index - 1];
            b_plus_1_imag[i] = alpha (B - index, l, r_squared[i]) * b_plus_1_imag[i] + b_plus_2_imag[i] + idata[B - l - index - 1];
            
            b_plus_2_real[i] = beta (B - index, l) * b_temp_real[i];
            b_plus_2_imag[i] = beta (B - index, l) * b_temp_imag[i];
         }       
      }
      else
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            b_plus_1_real[i] = rdata[B - l - index - 1];
            b_plus_1_imag[i] = idata[B - l - index - 1];
         }    
      }
   }
   
   long double* R_l_plus_1_l_vec = new long double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      R_l_plus_1_l_vec[i] = R_l_plus_1_l (l, r[i]);
   
   if ( B - l == 1 )
   {
      for ( int i = 0; i < 2 * B; ++i )
      {
         data_trans_real[i] = rdata[0] * R_l_plus_1_l_vec[i];
         data_trans_imag[i] = idata[0] * R_l_plus_1_l_vec[i];  
      }
      
      delete[] R_l_plus_1_l_vec;
      
      return;
   }
   
   long double* R_l_plus_2_l_vec = new long double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      R_l_plus_2_l_vec[i] = R_l_plus_2_l (l, r[i]);
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_trans_real[i] = (rdata[0] + b_plus_2_real[i]) * R_l_plus_1_l_vec[i] + b_plus_1_real[i] * R_l_plus_2_l_vec[i];
      data_trans_imag[i] = (idata[0] + b_plus_2_imag[i]) * R_l_plus_1_l_vec[i] + b_plus_1_imag[i] * R_l_plus_2_l_vec[i];
   }  
   
   delete[] R_l_plus_1_l_vec;
   delete[] R_l_plus_2_l_vec;
   
   if ( B - l - 2 >= 1 )
   {
      delete[] b_temp_real;
      delete[] b_temp_imag;
   }
     
   delete[] b_plus_1_real;
   delete[] b_plus_1_imag;
   
   delete[] b_plus_2_real;
   delete[] b_plus_2_imag;
   
   delete[] r_squared;
      
   return;
}
