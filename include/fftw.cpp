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

#include <math.h>

#include "fftw.h"
#include "fftw3.h"

using namespace std;

void fftw::dct_mdfd (int N, double* a, double* b)
{
   fftw_plan plan1 = fftw_plan_r2r_1d(N, a, a, FFTW_REDFT10, FFTW_ESTIMATE);
   
   fftw_execute(plan1);
   
   fftw_destroy_plan(plan1);
   
   fftw_plan plan2 = fftw_plan_r2r_1d(N, b, b, FFTW_REDFT10, FFTW_ESTIMATE);
   
   fftw_execute(plan2);
   
   fftw_destroy_plan(plan2);   
   
   a[0] /= 2.0 * N;
   b[0] /= 2.0 * N;
         
   for ( int i = 1; i < N; ++i )
   {
      a[i] /= N;
      b[i] /= N;
   }   
   
   return;
}

void fftw::dct_mdfd_adjoint (int N, double* a, double* b)
{
   for ( int i = 0; i < N; ++i )
   {
      a[i] /= N;
      b[i] /= N;
   }   
   
   fftw_plan plan1 = fftw_plan_r2r_1d(N, a, a, FFTW_REDFT01, FFTW_ESTIMATE);
   
   fftw_execute(plan1);
   
   fftw_destroy_plan(plan1);
   
   fftw_plan plan2 = fftw_plan_r2r_1d(N, b, b, FFTW_REDFT01, FFTW_ESTIMATE);
   
   fftw_execute(plan2);
   
   fftw_destroy_plan(plan2);   
   
   return;
}