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

#ifndef fsglft_h
#define fsglft_h

#include <string.h>
#include <fstream>

void fsglft (int B, double*** f_real, double*** f_imag, double*** f_hat_real, double*** f_hat_imag, std::string path);

void ifsglft (int B, double*** f_real, double*** f_imag, double*** f_hat_real, double*** f_hat_imag, std::string path);

long double R (int n, int l, long double r);

long double R_l_plus_1_l (int l, long double r);

long double R_l_plus_2_l (int l, long double r);

long double alpha (int n, int l, long double r_squared);

long double beta (int n, int l);

void DRT_Clenshaw (int B, int l, long double* r, long double* data_real, long double* data_imag, long double* data_trans_real, long double* data_trans_imag);

void iDRT_Clenshaw (int B, int l, long double* r, long double* data_real, long double* data_imag, long double* data_trans_real, long double* data_trans_imag);

#endif
