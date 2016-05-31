
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

function R_nl = R (n, l, r)

   r_squared = r .^ 2;
   
   temp = zeros (length (r), n - l);

   temp(:, 1) = 1.0;

   if ( n - l - 1 ~= 0 )
   
      temp(:, 2) = 1.5 + l - r_squared(:);

      for i = 2 : n - l - 1
       
        temp(:, i + 1) = (( 2 * i - 0.5 + l - r_squared(:) ) .* temp(:, i)     ...
                       + (- i + 0.5 - l) * temp(:, i - 1)) ...
                       / i;
   
      end
   
   end
   
   R_nl = temp(:, end);
   
   R_nl = R_nl .* r .^ l;
   
   factor = 2.0 ^ (n + 1.0) / 1.77245385090551588191; % sqrt (pi)
   
   for index = 2 * n - 1 : - 2 : 2
       
       factor = factor / index;
       
   end
   
   factor = factor * factorial (n - l - 1);
       
   R_nl = R_nl .* sqrt (factor);
   
return;
