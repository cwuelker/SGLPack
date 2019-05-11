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

#include "flt.h"

extern "C"
{
   #include <SpharmonicKit27/cospmls.h>
   #include <SpharmonicKit27/newFCT.h>
   #include <SpharmonicKit27/weights.h>
   #include <SpharmonicKit27/oddweights.h>   
}

void FLT::flt_seminaive (int bw, int m, double* data, double* data_transformed)
{
   const double pi = 3.1415926535897932384626433832795028841971;
   const double sqrt_pi = 1.77245385090551602729816748334;
   const double sqrt_2 = 1.41421356237309504880168872421;
   const double sqrt_2_pi = 2.50662827463100050241576528481;

   bool sign = false;

   if ( m < 0 )
   {
      sign = true;
      m = - m;
   }
   
   double* result = data_transformed;
   double* workspace = new double[8 * bw];
   double* cos_even = new double[bw];
   
   int tablesize = TableSize (m, bw);

   double* cos_pml_table = new double[tablesize];
   double* workspace_2 = new double[16 * bw];

   CosPmlTableGen (bw, m, cos_pml_table, workspace_2);

   // ------------------------------------------------------------------
   
   //for ( int index = 0; index < 2 * bw; ++index )
   //  data[index] /= sin ((2 * index + 1) * pi / (4.0 * bw));
   
   int i, j, n;
   const double *weights;
   double *weighted_data;
   double *cos_data;
   double *scratchpad;

   double *cos_odd;

   double eresult0, eresult1, eresult2, eresult3;
   double oresult0, oresult1;  
   double *pml_ptr_even, *pml_ptr_odd;

   int ctr;
   double d_bw;

   n = 2*bw;
   d_bw = (double) bw;

   cos_odd = cos_even + (bw/2);

   /* assign workspace */
   weighted_data = workspace;
   cos_data = workspace + (2*bw);
   scratchpad = cos_data + (2*bw); /* scratchpad needs (4*bw) */
   /* total workspace = (8 * bw) */
 
   /*
   need to apply quadrature weights to the data and compute
   the cosine transform: if the order is even, get the regular
   weights; if the order is odd, get the oddweights (those
   where I've already multiplied the sin factor in
   */

   if( (m % 2) == 0)
   {
      //weights = get_weights(bw);
   }
   else
   {
      //weights = get_oddweights(bw);
   }
 
   /*
   smooth the weighted signal
   */
   for ( i = 0; i < n    ; ++i )
     weighted_data[i] = data[ i ];// * weights[ i ];
  
   kFCT(weighted_data, cos_data, scratchpad, n, bw, 1);
  

   /* normalize the data for this problem */
   if (m == 0)
   {
     for (i=0; i<bw; i++)
       cos_data[i] *= (d_bw * 0.5);
   }
   else
   {
     for (i=0; i<bw; i++)
       cos_data[i] *= d_bw ;
   }

   cos_data[0] *= 2.0;


   /*** separate cos_data into even and odd indexed arrays ***/
   ctr = 0;
   for(i = 0; i < bw; i += 2)
   {
      cos_even[ctr] = cos_data[i];
      cos_odd[ctr]  = cos_data[i+1];
      //fprintf(stderr,"cos_even = %f\t cos_odd = %f\n", cos_data[i], cos_data[i+1]);
      ctr++;
   }
  
   /*
   do the projections; Note that the cos_pml_table has
   had all the zeroes stripped out so the indexing is
   complicated somewhat
   */

   /******** this is the original routine 
   for (i=m; i<bw; i++)
   {
   pml_ptr = cos_pml_table + NewTableOffset(m,i);

   if ((m % 2) == 0)
   {
   for (j=0; j<(i/2)+1; j++)
   result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
   }
   else
   {
   if (((i-m) % 2) == 0)
   {
   for (j=0; j<(i/2)+1; j++)
   result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
   }
   else
   {
   for (j=0; j<(i/2); j++)
   result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
   }
   } 
      
   toggle = (toggle+1) % 2;
   }

   end of original routine ***/

   /*** new version of the above for-loop; I'll
   try to remove some conditionals and move others;
   this may make things run faster ... I hope ! ***/

   /*
   this loop is unrolled to compute two coefficients
   at a time (for efficiency)
   */

   for(i = 0; i < (bw - m) - ( (bw - m) % 2 ); i += 2)
   {
      /* get ptrs to precomputed data */
      pml_ptr_even = cos_pml_table+ NewTableOffset(m, m + i);
      pml_ptr_odd  = cos_pml_table+ NewTableOffset(m, m + i + 1);

      eresult0 = 0.0; eresult1 = 0.0;
      oresult0 = 0.0; oresult1 = 0.0;

      for(j = 0; j < (((m + i)/2) + 1) % 2; ++j)
        {
          eresult0 += cos_even[j] * pml_ptr_even[j];
          oresult0 += cos_odd[j] * pml_ptr_odd[j];
        }
      for( ;  j < ((m + i)/2) + 1; j += 2)
        {
          eresult0 += cos_even[j] * pml_ptr_even[j];
          eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
          oresult0 += cos_odd[j]* pml_ptr_odd[j];
          oresult1 += cos_odd[j + 1]* pml_ptr_odd[j + 1];
        }

      /* save the result */
      result[i] =eresult0 + eresult1;
      result[i + 1] =oresult0 + oresult1;
   }

   /* if there are an odd number of coefficients to compute
   at this order, get that last coefficient! */
   if( ( (bw - m) % 2 ) == 1)
   {
      pml_ptr_even = cos_pml_table + NewTableOffset(m, bw - 1);

      eresult0 = 0.0; eresult1 = 0.0;
      eresult2 = 0.0; eresult3 = 0.0;

      for(j = 0; j < (((bw - 1)/2) + 1) % 4; ++j)
        {
          eresult0 += cos_even[j] * pml_ptr_even[j];
        }

      for( ;  j < ((bw - 1)/2) + 1; j += 4)
        {
          eresult0 += cos_even[j] * pml_ptr_even[j];
          eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
          eresult2 += cos_even[j + 2] * pml_ptr_even[j + 2];
          eresult3 += cos_even[j + 3] * pml_ptr_even[j + 3];
        }
      result[bw - m - 1] = eresult0 + eresult1 + eresult2+ eresult3;
   }
 
   // ------------------------------------------------------------------
    
   for ( int index = 0; index < bw - m; ++index )
   {
      if ( m == 0 )
      {
         result[index] /= sqrt_pi;
      }
      else
      {
         result[index] /= sqrt_2_pi;         
      }
   }
      
   if ( sign && m % 2 )
      for ( int index = 0; index < bw - m; ++index ) result[index] = - result[index];   
         
   delete[] workspace;
   delete[] workspace_2;
   delete[] cos_even;
   delete[] cos_pml_table;
   
   return;
}

void FLT::flt_seminaive_adjoint (int bw, int m, double* data, double* data_transformed)
{
   const double sqrt_pi = 1.77245385090551602729816748334;
   const double sqrt_2 = 1.41421356237309504880168872421;
   const double sqrt_2_pi = 2.50662827463100050241576528481;

   bool sign = false;

   if ( m < 0 )
   {
      sign = true;
      m = - m;
   }

   double* assoc_legendre_series = data;
   double* result = data_transformed;
   double* workspace = new double[5 * bw];

   double* sin_values = new double[2 * bw];

   for ( int i = 0; i < 2 * bw; ++i ) sin_values[i] = 1.0;

   int tablesize = TableSize (m, bw);

   double* tablespace = new double[tablesize];
   double* workspace_2 = new double[16 * bw];
   double* trans_cos_pml_table = new double[tablesize];

   CosPmlTableGen (bw, m, tablespace, workspace_2);
   Transpose_CosPmlTableGen(bw, m, tablespace, trans_cos_pml_table);

   // ------------------------------------------------------------------

   double *fcos; /* buffer for cosine series rep of result */
   double *trans_tableptr;
   double *scratchpad;
   double *assoc_offset;
   int i, j, n, rowsize;

   double *p;
   double fcos0, fcos1, fcos2, fcos3;

   fcos = workspace; /* needs bw */
   scratchpad = fcos + bw; /* needs (4*bw) */

   /* total workspace is 5*bw */

   trans_tableptr = trans_cos_pml_table;
   p = trans_cos_pml_table;

   /* main loop - compute each value of fcos

    Note that all zeroes have been stripped out of the
    trans_cos_pml_table, so indexing is somewhat complicated.
    */

   /***** original version of the loop (easier to see
    how things are being done)

    for (i=0; i<bw; i++)
    {
    if (i == (bw-1))
    {
    if ((m % 2) == 1)
    {
    fcos[bw-1] = 0.0;
    break;
    }
    }

    fcos[i] = 0.0;
    rowsize = Transpose_RowSize(i, m, bw);

    if (i > m)
    assoc_offset = assoc_legendre_series + (i - m) + (m % 2);
    else
    assoc_offset = assoc_legendre_series + (i % 2);

    for (j=0; j<rowsize; j++)
    fcos[i] += assoc_offset[2*j] * trans_tableptr[j];

    trans_tableptr += rowsize;
    }

    end of the original version ******/

   /*** this is a more efficient version of the above, original
    loop ***/

   for (i=0; i<=m; i++)
      {
      rowsize = Transpose_RowSize(i, m, bw);

      /* need to point to correct starting value of Legendre series */
      assoc_offset = assoc_legendre_series + (i % 2);

      fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

      for (j = 0; j < rowsize % 4; ++j)
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];

      for ( ; j < rowsize; j += 4){
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];
         fcos1 += assoc_offset[2*(j+1)] * trans_tableptr[j+1];
         fcos2 += assoc_offset[2*(j+2)] * trans_tableptr[j+2];
         fcos3 += assoc_offset[2*(j+3)] * trans_tableptr[j+3];
      }

      fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

      trans_tableptr += rowsize;
      }

   for (i=m+1; i<bw-1; i++)
   {

      rowsize = Transpose_RowSize(i, m, bw);

      /* need to point to correct starting value of Legendre series */
      assoc_offset = assoc_legendre_series + (i - m) + (m % 2);

      fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

      for (j = 0; j < rowsize % 4; ++j)
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];

      for ( ; j < rowsize; j += 4){
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];
         fcos1 += assoc_offset[2*(j+1)] * trans_tableptr[j+1];
         fcos2 += assoc_offset[2*(j+2)] * trans_tableptr[j+2];
         fcos3 += assoc_offset[2*(j+3)] * trans_tableptr[j+3];
      }
      fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

      trans_tableptr += rowsize;
   }

   if((m % 2) == 1)
      /* if m odd, no need to do last row - all zeroes */
      fcos[bw-1] = 0.0;
   else
      {
      rowsize = Transpose_RowSize(bw-1, m, bw);

      /* need to point to correct starting value of Legendre series */
      assoc_offset = assoc_legendre_series + (bw - 1 - m) + (m % 2);

      fcos0 = 0.0; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

      for(j = 0; j < rowsize % 4; ++j)
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];

      for ( ; j < rowsize; j += 4){
         fcos0 += assoc_offset[2*j] * trans_tableptr[j];
         fcos1 += assoc_offset[2*(j+1)] * trans_tableptr[j+1];
         fcos2 += assoc_offset[2*(j+2)] * trans_tableptr[j+2];
         fcos3 += assoc_offset[2*(j+3)] * trans_tableptr[j+3];
      }

      fcos[bw - 1] = fcos0 + fcos1 + fcos2 + fcos3;

      trans_tableptr += rowsize;
      }


   /* now we have the cosine series for the result,
    so now evaluate the cosine series at 2*bw Chebyshev nodes
    */

   ExpIFCT(fcos, result, scratchpad, 2*bw, bw, 1);

   /* if m is odd, then need to multiply by sin(x) at Chebyshev nodes */

   if ((m % 2) == 1)
      {
      n = 2*bw;

      for (j=0; j<n; j++)
         result[j] *= sin_values[j];
      }

   trans_tableptr = p;

   /* amscray */

   // ------------------------------------------------------------------

   if ( m == 0 )
      for (int i = 0; i < 2 * bw; ++i) 
         data_transformed[i] /= sqrt_2;
   
   if ( sign && m % 2 )
      for (int i = 0; i < 2 * bw; ++i) 
         data_transformed[i] = - data_transformed[i];

   for (int i = 0; i < 2 * bw; ++i) 
      data_transformed[i] /= sqrt_2_pi;

   delete[] workspace;
   delete[] workspace_2;
   delete[] tablespace;
   delete[] trans_cos_pml_table;
   delete[] sin_values;
   
   return;
}