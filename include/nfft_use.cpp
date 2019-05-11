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

#include <iostream>
#define NFFT_PRECISION_DOUBLE

#include "nfft_use.h"
#include "nfft3mp.h"

using namespace std;

void nfft_use::nfft (int N[3], int M, double*** data_real, double*** data_imag, double** points, complex<double>* result, int sigma, int q)
{
   int N_sigma[3];
   
   N_sigma[0] = sigma * N[0];
   N_sigma[1] = sigma * N[1];
   N_sigma[2] = sigma * N[2];
   
   NFFT (plan) p;
   NFFT (init_guru)(&p, 3, N, M, N_sigma, q, PRE_PHI_HUT | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

   int index = 0;
   
   for ( int i = 0; i < M; ++i )
   {
      for ( int j = 0; j < 3; ++j )
      {
         p.x[index] = points[i][j];
         
         ++index;
      }
   }

   index = 0;
   
   for (int i = 0; i < N[0]; ++i)
   {
      for (int j = 0; j < N[1]; ++j)
      {
         for (int s = 0; s < N[2]; ++s)
         {
            p.f_hat[index][0] = data_real[i][j][s];
            p.f_hat[index][1] = data_imag[i][j][s];

            ++index;
         }
      }
   }

   NFFT (precompute_one_psi)(&p);
   
   /** check for valid parameters before calling any trafo/adjoint method */
   const char* error_str = NFFT(check)(&p);
  
   if (error_str != 0)
      printf("Error in nfft module: %s\n", error_str);
   
   NFFT (trafo)(&p);

   for ( int i = 0; i < M; ++i )
   {
      result[i].real (p.f[i][0]);
      result[i].imag (p.f[i][1]);
   }

   NFFT (finalize)(&p);

   return;
}

void nfft_use::nfft_adjoint (int N[3], int M, double*** data_real, double*** data_imag, double** points, complex<double>* input, int sigma, int q)
{
   int N_sigma[3];
   
   N_sigma[0] = sigma * N[0];
   N_sigma[1] = sigma * N[1];
   N_sigma[2] = sigma * N[2];
   
   NFFT (plan) p;
   NFFT (init_guru)(&p, 3, N, M, N_sigma, q, PRE_PHI_HUT | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

   int index = 0;
   
   for ( int i = 0; i < M; ++i )
   {
      for ( int j = 0; j < 3; ++j )
      {
         p.x[index] = points[i][j];
         
         ++index;
      }
   }

   index = 0;
   
   for ( int i = 0; i < M; ++i )
   {
      p.f[i][0] = input[i].real ();
      p.f[i][1] = input[i].imag ();
   }
   
   NFFT (precompute_one_psi)(&p);
   
   /** check for valid parameters before calling any trafo/adjoint method */
   const char* error_str = NFFT(check)(&p);
  
   if (error_str != 0)
      printf("Error in nfft module: %s\n", error_str);
   
   NFFT (adjoint_3d)(&p);

   for (int i = 0; i < N[0]; ++i)
   {
      for (int j = 0; j < N[1]; ++j)
      {
         for (int s = 0; s < N[2]; ++s)
         {
            data_real[i][j][s] = p.f_hat[index][0];
            data_imag[i][j][s] = p.f_hat[index][1];

            ++index;
         }
      }
   }

   NFFT (finalize)(&p);

   return;
}