
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

function F = iDSGLFT (B, F_Fourier, path)

filename_cos_theta = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_cos_theta.txt', path, B, B);
filename_phi = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_phi.txt', path, B, B);
filename_N = sprintf ('%s/FSFT/precomp/%d/%d_100_digits_N.txt', path, B, B);

file_id_cos_theta = fopen (filename_cos_theta, 'r');
file_id_phi = fopen (filename_phi, 'r');
file_id_N = fopen (filename_N, 'r');

cos_theta = fscanf (file_id_cos_theta, '%f\n', 2 * B);
phi = fscanf (file_id_phi, '%f\n', 2 * B);
N = fscanf (file_id_N, '%f\n', B ^ 2);

fclose (file_id_cos_theta);
fclose (file_id_phi);
fclose (file_id_N);

r = half_range_Hermite_quadrature_abcissae_precomp (2 * B, path);

F = zeros (2 * B, 2 * B, 2 * B);

for n = 1 : B
     
 for i = 0 : 2 * B - 1
  for k = 0 : 2 * B - 1        
    
   counter = 0;     
      
   for l = 0 : n - 1
    
    P_l_all_m_all_theta = P_all_m_all_t (l, cos_theta);   
    
     for m = - l : l
     
      counter = counter + 1;  
      
      for j = 0 : 2 * B - 1
     
          R_nl_i = R_mex (n, l, r(i + 1));           
           
          F(j + 1, k + 1, i + 1) = F(j + 1, k + 1, i + 1) + F_Fourier(n, l + 1, m + B) * R_nl_i * N(counter) * P_l_all_m_all_theta(l + m + 1, j + 1) * exp(1j * m * phi(k + 1));
           
       end
      end
     end      
      
  end
 end
end

return;


function P_l_all_m_all_t = P_all_m_all_t (l, t_vec)

   P_l_all_m_all_t = zeros (2 * l + 1, length (t_vec));

   P = legendre (l, t_vec);
   
   for m = - l : l
   
      if ( m < 0 )
       
         P_l_all_m_all_t(m + l + 1, :) = (-1.0) ^ m * factorial (l + m) * P(1 - m, :) / factorial (l - m);
       
      else

         P_l_all_m_all_t(m + l + 1, :) = P(1 + m, :);
       
      end

   end
   
return;
