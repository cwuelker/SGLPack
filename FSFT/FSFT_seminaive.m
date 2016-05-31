
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

function F_Fourier = FSFT_seminaive (B, P, weights, N, F)

F_Fourier = zeros (2 * B - 1, B);

for j = 0 : 2 * B - 1
    
   F(j + 1, :) = fft (F(j + 1, :));
    
end

for k = 0 : 2 * B - 1
    
   F(:, k + 1) = dct (weights .* F(:, k + 1));
   
end

l_index = 1;

for m = - B + 1 : B - 1
 
 for l = abs (m) : B - 1
     
    start_index = 2 - mod (l, 2); 
       
    if ( m < 0 )
            
     if ( mod (m, 2) )   
       
        F_Fourier(m + B, l + 1) = P(l_index, start_index : 2 : end) * F(start_index : 2 : end, 2 * B + m + 1);
     
     else
         
        F_Fourier(m + B, l + 1) = P(l_index, 1 : l + 1) * F(1 : l + 1, 2 * B + m + 1); 
         
     end
       
    else
          
     if ( mod (m, 2) )   
        
        F_Fourier(m + B, l + 1) = P(l_index, start_index : 2 : end) * F(start_index : 2 : end, m + 1);
     
     else
         
        F_Fourier(m + B, l + 1) = P(l_index, 1 : l + 1) * F(1 : l + 1, m + 1);
         
     end        
        
    end
    
    l_index = l_index + 1;
            
 end
   
end

counter = 0;

for l = 0 : B - 1
 for m = - l : l
   
   counter = counter + 1;  
     
   F_Fourier(m + B, l + 1) = pi * F_Fourier(m + B, l + 1) * N(counter) / B;
    
 end
end

return;
