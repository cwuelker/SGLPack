
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

function [abcissae, weights] = half_range_Hermite_quadrature (n)

format long;

n = n - 1;

gamma_vec = zeros (n + 2, 1);
p_zero_vec = zeros (n + 2, 1);

alpha_vec = zeros (n + 1, 1);
beta_vec = zeros (n + 1, 1);

alpha_vec (1) = - 1.0 / sqrt (pi);

gamma_vec (1) = 0.50 * sqrt (pi);
gamma_vec (2) = 0.25 * sqrt (pi) - 0.5 / sqrt (pi);

p_zero_vec (1) = 1.0;
p_zero_vec (2) = - 1.0 / sqrt (pi);

for index = 3 : n + 2

   alpha_vec(index - 1) = - (0.5 / gamma_vec(index - 1)) * p_zero_vec(index - 1)^2;
   beta_vec(index - 1) = - gamma_vec(index - 1) / gamma_vec(index - 2);
   
   p_zero_vec(index) = alpha_vec(index - 1) * p_zero_vec(index - 1) + beta_vec(index - 1) * p_zero_vec(index - 2);
   
   gamma_vec(index) = 0.5 * (index - 1) * gamma_vec(index - 1) + 0.5 * p_zero_vec(index) * p_zero_vec(index - 1);

end

J_n = zeros (n + 1);

for index = 1 : n + 1
    
   J_n(index, index) = - alpha_vec(index);
    
end

for index = 2 : n + 1
    
   J_n(index, index - 1) = sqrt (- beta_vec(index));
   J_n(index - 1, index) = J_n(index, index - 1); 
   
end

abcissae = sort (eig(J_n));

% -------------------------------------------------------------------------

weights = zeros (n + 1, 1);

[value_vec, derivative_vec] = precomp (abcissae, alpha_vec, beta_vec);

for index = 1 : n + 1
    
   weights(index) = gamma_vec(n + 1) / (derivative_vec (index) * value_vec(index));
    
end

return;


function [value_vec, derivative_vec] = precomp (x_vec, alpha_vec, beta_vec)

   n = length (alpha_vec);
   
   p_k_vec = zeros (1, n);
   p_k_minus_1_vec = zeros (1, n);

   p_k_vec(1) = - 1.0 / sqrt (pi);
   p_k_vec(2) = 1.0;
   p_k_minus_1_vec(1) = 1.0;
   
   pd_k_vec = zeros (1, n);
   pd_k_minus_1_vec = zeros (1, n);
   
   pd_k_vec(1) = 1.0;
   
   for index = 1 : n - 1
   
      pd_k_vec_temp = pd_k_vec;
    
      pd_k_vec = p_k_vec + [0.0, pd_k_vec(1 : end - 1)] + alpha_vec(index + 1) * pd_k_vec + beta_vec(index + 1) * pd_k_minus_1_vec;
      pd_k_minus_1_vec = pd_k_vec_temp;
      
      if ( index ~= n - 1 )
      
         p_k_vec_temp = p_k_vec;
    
         p_k_vec = [0.0, p_k_vec(1 : end - 1)] + alpha_vec(index + 1) * p_k_vec + beta_vec(index + 1) * p_k_minus_1_vec;
         p_k_minus_1_vec = p_k_vec_temp;
         
      end
   
   end
   
   value_vec = zeros (1, n);
   derivative_vec = zeros (1, n);
   
   for index_1 = 1 : n
    for index_2 = 1 : n
    
     value_vec(index_1) = value_vec(index_1) + p_k_vec(index_2) * x_vec(index_1) ^ (index_2 - 1);   
     derivative_vec(index_1) = derivative_vec(index_1) + pd_k_vec(index_2) * x_vec(index_1) ^ (index_2 - 1); 
 
    end
   end
   
return;
