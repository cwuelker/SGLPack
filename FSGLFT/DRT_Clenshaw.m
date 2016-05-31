
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

function data_trans = DRT_Clenshaw (B, l, r, data)

data_trans = zeros (B - l, 1);
data_temp = zeros (4 * B, 1);

for i = 0 : 2 * B - 1
    
 data_temp(2 * B + i + 1) = R_l_plus_1_l (l, r(i + 1));  
    
end

data_trans(1) = data_temp(2 * B + 1 : end)' * data;

if ( B - l == 1 )
    
 return;
 
end

for i = 0 : 2 * B - 1
    
 data_temp(i + 1) = R_l_plus_2_l (l, r(i + 1));  
    
end

data_temp(1 : 2 * B) = data .* data_temp(1 : 2 * B);

data_trans(2) = sum (data_temp(1 : 2 * B));

if ( B - l == 2 )
    
 return;
 
end

data_temp(2 * B + 1 : end) = beta (l + 2, l) * data .* data_temp(2 * B + 1 : end);

for index = 1 : B - l - 2
    
 if ( index < B - l - 2 )   
 
  temp = data_temp(1 : 2 * B);   
  
 end 
 
 data_temp(1 : 2 * B) = alpha (l + 1 + index, l, r_squared) .* data_temp(1 : 2 * B) + data_temp(2 * B + 1 : end);
 
 data_trans(2 + index) = sum (data_temp(1 : 2 * B));
 
 if ( index < B - l - 2 )
     
  data_temp(2 * B + 1 : end) = beta (l + 2 + index, l) * temp;
 
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


function R_val = R_l_plus_1_l (l, r)

 R_val = R_mex (l + 1, l, r) * exp (- r ^ 2);
 
return;


function R_val = R_l_plus_2_l (l, r)

 R_val = R_mex (l + 2, l, r) * exp (- r ^ 2);

return;


function alpha_val = alpha (n, l, r_squared)

 alpha_val = (2 * n - l - 0.5 - r_squared) / sqrt ((n + 0.5) * (n - l));

return;


function beta_val = beta (n, l)

 beta_val = - sqrt ((n - 0.5) * (n - l - 1) / ((n + 0.5) * (n - l)));

return;
