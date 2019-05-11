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

#include <fsft.h>

extern "C"
{
   #include <FST_semi_memo.h>
   #include <cospmls.h>
}
 
void fsft (int B, double** f_real, double** f_imag, double** f_hat_real, double** f_hat_imag)
{   
   const int spharmonic_tablesize = Reduced_Naive_TableSize (B, 0) + Reduced_SpharmonicTableSize (B, 0);
   
   double* resultspace = new double[spharmonic_tablesize];
   
   double* workspace = new double[16 * B];
   double* workspace_2 = new double[(8 * B * B) + (29 * B)];
   
   double* f_hat_real_fsft = new double[B * B];
   double* f_hat_imag_fsft = new double[B * B];
   
   double** seminaive_pml_table = SemiNaive_Naive_Pml_Table (B, B, resultspace, workspace); 
   
   double* f_real_fsft = new double[4 * B * B];
   double* f_imag_fsft = new double[4 * B * B];
   
   int* perm = new int[2 * B];
   
   for ( int index = 0; index <= B; ++index )
      perm[index] = B - index;
   
   for ( int index = 0; index < B - 1; ++index )
      perm[index + B + 1] = 2 * B - 1 - index;
   
   for ( int j = 0; j < 2 * B; ++j )
   {
      for ( int k = 0; k < 2 * B; ++k )
      {
         f_real_fsft[2 * B * j + perm[k]] = f_real[j][k];
         f_imag_fsft[2 * B * j + perm[k]] = f_imag[j][k];
      }
   }
   
   FST_semi_memo (f_real_fsft, f_imag_fsft, f_hat_real_fsft, f_hat_imag_fsft, 2 * B, seminaive_pml_table, workspace_2, 0, B); 
   
   for ( int l = 0; l < B; ++l )
   {       
      for ( int m = - l; m <= l; ++m )
      {        
         if ( m < 0 )
         {
            int offset = 0;
         
            for ( unsigned index = 0; index < - m; ++index )
               offset += B - index;

            f_hat_real[l][m + B - 1] = f_hat_real_fsft[l + m + offset];
            f_hat_imag[l][m + B - 1] = f_hat_imag_fsft[l + m + offset];
         }
         else if ( m == 0 )
         {
            f_hat_real[l][B - 1] = f_hat_real_fsft[l];
            f_hat_imag[l][B - 1] = f_hat_imag_fsft[l];
         }
         else  
         {
            int offset = 0;
            
            for ( unsigned index = 1; index < m; ++index )
               offset += B - index;
            
            f_hat_real[l][m + B - 1] = f_hat_real_fsft[B * (B - 1) + l - offset];
            f_hat_imag[l][m + B - 1] = f_hat_imag_fsft[B * (B - 1) + l - offset];
         }
      }
   }
   
   delete[] perm;
   
   delete[] f_real_fsft;
   delete[] f_imag_fsft;
   
   delete[] f_hat_real_fsft;
   delete[] f_hat_imag_fsft;
   
   delete[] resultspace;
   
   delete[] workspace;
   delete[] workspace_2;
   
   delete[] seminaive_pml_table;

   return;
}

void ifsft (int B, double** f_real, double** f_imag, double** f_hat_real, double** f_hat_imag)
{   
   const int spharmonic_tablesize = Reduced_Naive_TableSize (B, 0) + Reduced_SpharmonicTableSize (B, 0);
   
   double* resultspace = new double[spharmonic_tablesize];
   double* resultspace_2 = new double[spharmonic_tablesize];
   
   double* workspace = new double[16 * B];
   double* workspace_2 = new double[8 * B * B + 32 * B];
   
   double* f_hat_real_fsft = new double[B * B];
   double* f_hat_imag_fsft = new double[B * B];
      
   double** seminaive_pml_table = SemiNaive_Naive_Pml_Table (B, B, resultspace, workspace); 
   double** seminaive_pml_table_inv = Transpose_SemiNaive_Naive_Pml_Table (seminaive_pml_table, B, B, resultspace_2, workspace);
      
   for ( int l = 0; l < B; ++l )
   {       
      for ( int m = - l; m <= l; ++m )
      {        
         if ( m < 0 )
         {
            int offset = 0;
         
            for ( unsigned index = 0; index < - m; ++index )
               offset += B - index;

            f_hat_real_fsft[l + m + offset] = f_hat_real[l][m + l];
            f_hat_imag_fsft[l + m + offset] = f_hat_imag[l][m + l];
         }
         else if ( m == 0 )
         {
            f_hat_real_fsft[l] = f_hat_real[l][l];
            f_hat_imag_fsft[l] = f_hat_imag[l][l];
         }
         else  
         {
            int offset = 0;
            
            for ( unsigned index = 1; index < m; ++index )
               offset += B - index;
            
            f_hat_real_fsft[B * (B - 1) + l - offset] = f_hat_real[l][m + l];
            f_hat_imag_fsft[B * (B - 1) + l - offset] = f_hat_imag[l][m + l];
         }
      }
   }
   
   double* f_real_fsft = new double[4 * B * B];
   double* f_imag_fsft = new double[4 * B * B];
   
   int* perm = new int[2 * B];
   
   for ( int index = 0; index <= B; ++index )
      perm[index] = B - index;
   
   for ( int index = 0; index < B - 1; ++index )
      perm[index + B + 1] = 2 * B - 1 - index;
   
   InvFST_semi_memo (f_hat_real_fsft, f_hat_imag_fsft, f_real_fsft, f_imag_fsft, 2 * B, seminaive_pml_table_inv, workspace_2, 0, B);
   
   for ( int j = 0; j < 2 * B; ++j )
   {
      for ( int k = 0; k < 2 * B; ++k )
      {
         f_real[j][k] = f_real_fsft[2 * B * j + perm[k]];
         f_imag[j][k] = f_imag_fsft[2 * B * j + perm[k]];
      }
   }   
   
   delete[] perm;
   
   delete[] f_real_fsft;
   delete[] f_imag_fsft;
   
   delete[] f_hat_real_fsft;
   delete[] f_hat_imag_fsft;
   
   delete[] resultspace;
   delete[] resultspace_2;
   
   delete[] workspace;
   delete[] workspace_2;
   
   delete[] seminaive_pml_table;
   delete[] seminaive_pml_table_inv;
   
   return;
}
