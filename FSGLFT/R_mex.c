//
//        SGLPack - Fast spherical Gauss-Laguerre Fourier transforms        
//   
//  
//   Contact: Christian Wuelker, M.Sc.
//            wuelker@math.uni-luebeck.de
//  
//   Copyright 2016  Christian Wuelker, Juergen Prestin
//
//
//     SGLPack is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 3 of the License, or
//     (at your option) any later version.
//  
//     SGLPack is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

#include "mex.h"
#include <math.h>

int factorial (int n)
{
   return (n == 1 || n == 0) ? 1 : factorial (n - 1) * n;
}

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   
   double* R, * r, * n, * l;
   
   R = mxGetPr (plhs[0]);
   
   n = mxGetPr (prhs[0]);
   l = mxGetPr (prhs[1]);
   
   r = mxGetPr (prhs[2]);
   
   double R_val;
   
   const double r_val = *r;
   
   const int n_val = (int)(*n);
   const int l_val = (int)(*l);
   
   const int n_l_difference = n_val - l_val;
   
   const double r_val_squared = r_val * r_val;
   
   double temp[n_l_difference];
   
   unsigned index;
   
   temp[0] = 1.0;
   
   if ( n_l_difference > 1 )
   {
      temp[1] = 1.5 + l_val - r_val_squared;
   
      for ( index = 2; index < n_l_difference; ++index )
      {
         temp[index] = 0.0;
      }
   
      for ( index = 2; index < n_l_difference; ++index )
      {
          temp[index] = (((2.0 * index - 0.5 + l_val - r_val_squared) * temp[index - 1]) + ((0.5 - index - l_val) * temp[index - 2])) / index;
      }
   }
   
   R_val = temp[n_l_difference - 1];
  
   R_val *= pow (r_val, l_val);
   
   double factor = pow (2.0, n_val + 1.0) / 1.77245385090551588191;
   
   for ( index = 2 * n_val - 1; index > 1; index -= 2 )
   {
      factor /= index; 
   }
   
   for ( index = 2; index < n_l_difference; ++index )
   {
      factor *= index; 
   }
  
   R_val *= sqrt (factor);
   
   R[0] = R_val;
   
   return;
}

