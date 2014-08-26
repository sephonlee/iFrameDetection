function [J , t_est] = prgdyn( matD , Kmax);
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>

%%%%%%%%%%%% Programmation dynamique
% calcule les contrastes par dimension jusqu'a Kmax
% et les lieux des ruptures pour chaque dimension 

[N , N] = size(matD);

I   = repmat(Inf,Kmax,N);
t   = zeros(Kmax-1,N);
I(1,:) = matD(1,:);

if Kmax > 2,
   for k = 2:Kmax-1,
      for L = k:N
          [I(k,L),t(k-1,L)] = min(I(k-1,1:L-1) + matD(2:L,L)');
      end
   end
   [I(Kmax,N) , t(Kmax-1,N)] = min(I(Kmax-1,1:N-1) + matD(2:N,N)');
end

if Kmax==2    
   [I(Kmax,N) , t(Kmax-1,N)] = min(I(Kmax-1,1:N-1)+matD(2:N,N)');
end

J = I(:,N);

% *** Compute the change-point instants ***
t_est = diag(N*ones(Kmax,1));

for K = 2:Kmax,
    for k = K-1:-1:1,
        t_est(K,k) = t(k,t_est(K,k+1));
    end
end
