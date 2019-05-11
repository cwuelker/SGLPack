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

#include "clenshaw.h"

#include <cmath>

using namespace std;

long double clenshaw::R (int n, int l, long double r)
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

long double clenshaw::R_l_plus_1_l (int l, long double r)
{
   return R (l + 1, l, r);
}

long double clenshaw::R_l_plus_2_l (int l, long double r)
{
   return R (l + 2, l, r);
}

long double clenshaw::alpha (int n, int l, long double r_squared)
{
   return (2 * n - l - 0.5 - r_squared) / (long double)sqrt ((n + 0.5) * (n - l));
}

long double clenshaw::beta (int n, int l)
{
   return - (long double)sqrt ((long double)(n - 0.5) * (long double)(n - l - 1) / (long double)((n + 0.5) * (n - l)));
}

void clenshaw::DRT_Clenshaw (int B, int l, long double* r, long double* rdata, long double* idata, long double* data_trans_real, long double* data_trans_imag)
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

void clenshaw::DRT_Clenshaw_adjoint (int B, int l, long double* r, long double* rdata, long double* idata, long double* data_trans_real, long double* data_trans_imag)
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