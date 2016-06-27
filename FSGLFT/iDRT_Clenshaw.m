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

function data_trans = iDRT_Clenshaw (B, l, r, data)

r_squared = r .^ 2;

b_vec_plus_1 = zeros (2 * B, 1);
b_vec_plus_2 = zeros (2 * B, 1);

for index = 0 : B - l - 2  
  
 if ( index > 0 )
 
  b_vec_temp = b_vec_plus_1;
 
 end
  
 if ( index > 0 )
     
  b_vec_plus_1 = alpha (B - index, l, r_squared) .* b_vec_plus_1 + b_vec_plus_2 + data(B - l - index);
  
 else
    
  b_vec_plus_1 = data(B - l - index) * ones (2 * B, 1);   
     
 end
 
 if ( index > 0 )
     
  b_vec_plus_2 = beta (B - index, l) * b_vec_temp;
 
 end
  
end

R_l_plus_1_l_vec = zeros (2 * B, 1);

for i = 0 : 2 * B - 1
   
 R_l_plus_1_l_vec(i + 1) = R_l_plus_1_l (l, r(i + 1));
    
end

if ( B - l == 1 )
   
 data_trans = data(1) * R_l_plus_1_l_vec;
    
 return;   
    
end    
    
R_l_plus_2_l_vec = zeros (2 * B, 1);

for i = 0 : 2 * B - 1
   
 R_l_plus_2_l_vec(i + 1) = R_l_plus_2_l (l, r(i + 1));
    
end

data_trans = (data(1) + b_vec_plus_2) .* R_l_plus_1_l_vec + b_vec_plus_1 .* R_l_plus_2_l_vec;

return;


function R_val = R_l_plus_1_l (l, r)

 R_val = R_mex (l + 1, l, r);
 
return;


function R_val = R_l_plus_2_l (l, r)

 R_val = R_mex (l + 2, l, r);

return;


function alpha_val = alpha (n, l, r_squared)

 alpha_val = (2 * n - l - 0.5 - r_squared) / sqrt ((n + 0.5) * (n - l));

return;


function beta_val = beta (n, l)

 beta_val = - sqrt ((n - 0.5) * (n - l - 1) / ((n + 0.5) * (n - l)));

return;
