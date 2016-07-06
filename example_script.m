%
%        SGLPack - Fast spherical Gauss-Laguerre Fourier transforms        
%   
%  
%   Contact: Christian Wuelker, M.Sc.
%            wuelker@math.uni-luebeck.de
%  
%   Copyright 2016  Christian Wuelker, Juergen Prestin
%
%
%     SGLPack is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
%  
%     SGLPack is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%  
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%  
%  
%   Commercial use is absolutely prohibited.
%  
%   See the accompanying LICENSE file for details.

% -------------------------------------------------------------------------
fprintf ('\nSGLPack - Fast spherical Gauss-Laguerre Fourier transforms\n');
fprintf ('(Copyright 2016 Christian Wuelker, Juergen Prestin)\n\n');

% path to SGLPack
path = '/home/wuelker/www/SGLPack-master';

% add half-range Hermite quadrature and fast spherical Fourier transforms
% to path
% -------------------------------------------------------------------------
half_range_Hermite_quadrature_path = sprintf ('%s/half-range_Hermite_quadrature', path);
addpath (half_range_Hermite_quadrature_path);

FSFT_path = sprintf ('%s/FSFT', path);
addpath (FSFT_path);

FSGLFT_path = sprintf ('%s/FSGLFT', path);
addpath (FSGLFT_path);
% -------------------------------------------------------------------------

% set bandwidth
B = 32;

fprintf ('>> Bandwidth B = %d\n', B);

% read in precomputed data for half-range Hermite quadrature
% -------------------------------------------------------------------------
Hermite_r = half_range_Hermite_quadrature_abcissae_precomp (2 * B, path);
Hermite_a = half_range_Hermite_quadrature_adapted_weights_precomp (2 * B, path);
% -------------------------------------------------------------------------

% read in precomputed data for fast spherical Fourier transforms
% -------------------------------------------------------------------------
FSFT_P = zeros (B ^ 2, 2 * B);

filename_FSFT_weights = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_weights.txt', path, B, B);
filename_FSFT_N = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_N.txt', path, B, B);
filename_FSFT_N_sort = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_N_sort.txt', path, B, B);
filename_FSFT_P = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_P_seminaive.txt', path, B, B);

file_id_FSFT_weights = fopen (filename_FSFT_weights, 'r');
file_id_FSFT_N = fopen (filename_FSFT_N, 'r');
file_id_FSFT_N_sort = fopen (filename_FSFT_N_sort, 'r');
file_id_FSFT_P = fopen (filename_FSFT_P, 'r');

FSFT_weights = fscanf (file_id_FSFT_weights, '%f\n', 2 * B);
FSFT_N = fscanf (file_id_FSFT_N, '%f\n', B ^ 2);
FSFT_N_sort = fscanf (file_id_FSFT_N_sort, '%f\n', B ^ 2);
FSFT_P(:) = fscanf (file_id_FSFT_P, '%e\n', 2 * B ^ 3);

fclose (file_id_FSFT_weights);
fclose (file_id_FSFT_N);
fclose (file_id_FSFT_N_sort);
fclose (file_id_FSFT_P);
% -------------------------------------------------------------------------

% generate random SGL Fourier coefficients
% -------------------------------------------------------------------------
fprintf ('<< Generating random SGL Fourier coefficients... ');

F_Fourier = zeros (B, B, 2 * B - 1);

for m = - B + 1 : B - 1
 for l = abs (m) : B - 1
  for n = l + 1 : B
     
    F_Fourier(n, l + 1, m + B) = unifrnd (-1.0, 1.0) + 1j * unifrnd (-1.0, 1.0);
    
  end
 end
end

fprintf ('done.\n');
% -------------------------------------------------------------------------

fprintf ('<< Performing inverse fast SGL Fourier transform (iFSGLFT) to reconstruct function values from SGL Fourier coefficients... ');

tic;

% perform inverse fast SGL Fourier transform to reconstruct function values
F_recon = iFSGLFT_Clenshaw (B, FSFT_P', FSFT_N_sort, Hermite_r, F_Fourier);

runtime_inverse = toc;

fprintf ('done.\n');
fprintf ('<< Performing fast SGL Fourier transform (FSGLFT) to reconstruct SGL Fourier coefficients from reconstructed function values... ');

tic;

% perform fast SGL Fourier transform to reconstruct original SGL Fourier
% coefficients
F_Fourier_recon = FSGLFT_Clenshaw (B, FSFT_P, FSFT_weights, FSFT_N, Hermite_r, Hermite_a, F_recon);

runtime_forward = toc;

fprintf ('done.\n');

runtime = runtime_inverse + runtime_forward;

fprintf ('<< Total runtime: %f s (%d%% forward, %d%% inverse)\n', runtime, round (100.0 * runtime_forward / runtime), round (100.0 * runtime_inverse / runtime));

% determine maximum and relative transformation error
% -------------------------------------------------------------------------
max_error = 0.0;
rel_error = 0.0;

for m = - B + 1 : B - 1
 for l = abs (m) : B - 1
  for n = l + 1 : B
     
     if ( abs (F_Fourier(n, l + 1, m + B) - F_Fourier_recon(n, l + 1, m + B)) >= max_error )
         
        max_error = abs (F_Fourier(n, l + 1, m + B) - F_Fourier_recon(n, l + 1, m + B));
         
     end
     
     if ( (abs (F_Fourier(n, l + 1, m + B) - F_Fourier_recon(n, l + 1, m + B)) / abs (F_Fourier(n, l + 1, m + B))) >= rel_error )
         
        rel_error = abs (F_Fourier(n, l + 1, m + B) - F_Fourier_recon(n, l + 1, m + B)) / abs (F_Fourier(n, l + 1, m + B));
         
     end     
     
  end
 end
end

fprintf ('<< Maximum (relative) transformation error: %e (%e, respectively)\n\n', max_error, rel_error);
% -------------------------------------------------------------------------

% garbage collect
% -------------------------------------------------------------------------
clear file_id_FSFT_N file_id_FSFT_N_sort file_id_FSFT_P file_id_FSFT_weights
clear filename_FSFT_N filename_FSFT_N_sort filename_FSFT_P filename_FSFT_weights
clear FSFT_N FSFT_N_sort FSFT_P FSFT_path FSFT_weights FSGLFT_path
clear half_range_Hermite_quadrature_path Hermite_a Hermite_r n l m path
% -------------------------------------------------------------------------
