
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

function adapted_weights = half_range_Hermite_quadrature_adapted_weights_precomp (n, path)

filename_adapted_weights = sprintf ('%s/half-range_Hermite_quadrature/precomp/%d_100_digits_adapted_weights.txt', path, n);

file_id_adapted_weights = fopen (filename_adapted_weights, 'r');

adapted_weights = fscanf (file_id_adapted_weights, '%f\n', n);

fclose (file_id_adapted_weights);

return;
