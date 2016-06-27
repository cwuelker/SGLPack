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

function F = iFSFT_seminaive (B, P, N_sort, F_Fourier)

F = zeros (2 * B, 2 * B);

l_index = 1;

counter = 0;

for m = - B + 1 : B - 1
 
 for l = abs (m) : B - 1
     
    start_index = 2 - mod (l, 2); 
    
    counter = counter + 1;
    
    if ( m < 0 )
    
     if ( mod (m, 2) )  
        
        F(start_index : 2 : end, 2 * B + m + 1) = F(start_index : 2 : end, 2 * B + m + 1) + F_Fourier(m + B, l + 1) * N_sort(counter) * P(start_index : 2 : end, l_index);
    
     else
         
        F(1 : l + 1, 2 * B + m + 1) = F(1 : l + 1, 2 * B + m + 1) + F_Fourier(m + B, l + 1) * N_sort(counter) * P(1 : l + 1, l_index); 
         
     end
       
    else
       
     if ( mod (m, 2) )   
        
        F(start_index : 2 : end, m + 1) = F(start_index : 2 : end, m + 1) + F_Fourier(m + B, l + 1) * N_sort(counter) * P(start_index : 2 : end, l_index);
        
     else
         
        F(1 : l + 1, m + 1) = F(1 : l + 1, m + 1) + F_Fourier(m + B, l + 1) * N_sort(counter) * P(1 : l + 1, l_index);
         
     end
     
    end
     
    l_index = l_index + 1;
    
 end
 
end

for k = 0 : 2 * B - 1
   
   if ( k ~= B )

      F(:, k + 1) = idct (F(:, k + 1));
      
   end
      
end

for j = 0 : 2 * B - 1
    
   F(j + 1, :) = 2.0 * B * ifft (F(j + 1, :));
    
end

return;
