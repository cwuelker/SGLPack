
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

function F_Fourier = FSGLFT (B, FSFT_P, FSFT_weights, FSFT_N, Hermite_r, Hermite_a, F)

F_FSFT = zeros (2 * B, 2 * B - 1, B);
F_Fourier = zeros (B, B, 2 * B - 1);

for i = 0 : 2 * B - 1
    
   F_FSFT(i + 1, :, :) = FSFT_seminaive (B, FSFT_P, FSFT_weights, FSFT_N, F(:, :, i + 1));
    
end

for m = - B + 1 : B - 1
 for l = abs (m) : B - 1

    F_Fourier(l + 1 : end, l + 1, m + B) = R_mat (B, l, Hermite_r) * (Hermite_a .* F_FSFT(:, m + B, l + 1));

 end
end

return;


function R_mat_l = R_mat (B, l, r)

   R_mat_l = zeros (B - l, 2 * B);

   for n = l + 1 : B
    for i = 0 : 2 * B - 1
       
     R_mat_l(n - l, i + 1) = R_mex (n, l, r(i + 1)) * exp (- r(i + 1) ^ 2);

    end
   end

return;
