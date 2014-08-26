function [Y s]=DataGener(t,reg_rupt,reg_val,sig_rupt,sig_val,loi_err)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genere data with iid errors, for a regression function and a noise-level given as an input
%
% Input:
% t	: vector of observation instants
% reg_rupt : vector of the breakpoints of s: r_1 = 0 < r_2 < ... < r_K=1 
% reg_val : vecteur of the values of s on each of these intervals: the i-th value is the value of s over [r_i,r_{i+1}[, except for i=1 (value on (-\infty,r_1)) and i=K-1 (value on [r_{K-1},+\infty)) 
% sig_rupt : vector of the breakpoints of sigma (same convention as for s)
% sig_val : vecteur of the values of sigma on each of these intervals (same convention as for s)
% loi_err : integer for choosing the distribution of the errors
%		1 for Gaussian errors
%		2 for Exponential errors
% loi_err : nombre indiquant quelle est la loi des erreurs. 
%       1=Gaussien. 2=Poisson symetrise. 3=Laplace. 4=Exponentielle (asymetrique).
%
% Output:
% Y : vector of size n=numel(t) of the observations Y_1, ..., Y_n
% s : vector of size n of the values s(t_1), ..., s(t_n) of the regression function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=numel(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
switch loi_err
    case 1
        % Gaussian errors
        epsilon=randn(1,n);
    case 2
        % Exponentielle asymetrique (+2)
        tmp_exp=exprnd(1,1,n);
        epsilon=(tmp_exp-1)/1;
    otherwise
        epsilon=NaN*ones(1,n);
        fprintf('Unknown distribution for the errors in DataGener.m\n');
end

% Multiply by the values of the noise-level
Noise=eval_histo(sig_val,sig_rupt,t).*epsilon;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=eval_histo(reg_val,reg_rupt,t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=s+Noise;





