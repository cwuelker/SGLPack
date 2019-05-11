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

#ifndef nfsglft_h
#define nfsglft_h

void ndsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag);
void nfsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q);

void ndsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag);
void nfsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q);

void infsglft_golub_van_loan (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences);
void infsglft_bjoerck (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences);

long double R (int n, int l, long double r);

double M_P (int l, int m, double xi);
void Y (int l, int m, double theta, double phi, double* result);

#endif