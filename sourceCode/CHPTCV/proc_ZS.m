function [D_ZS, mu_ZS, rupt_ZS, crit2ZS, rupt_1ZS, mu_1ZS]=proc_ZS(Data, Dimmax, infos)
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
% Procedure [ZS]
%
% Input: 
% Data : vector of observations Y_1, ..., Y_n
% Dimmax : maximal dimension considered
% infos: structure with additional informations
%       infos.delta= (minimal size of a segment)-1
%           delta=0 : all models are kept
%           delta=1 : models with at least one segment of the form
%           [t_i,t_{i+1}) are removed
%           defaut: delta=1
%
% Output:
%
% D_ZS : dimension selected by the procedure [ZS] with sigma estimated ('ZS' in the paper)
% mu_ZS : values of the final estimator at t_1, ..., t_n
% rupt_ZS : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2ZS : vector of estimated values of the risk (with the ZS criterion) 
%           for all values of the dimension between 1 and Dimmax 
%
% rupt_1ZS : matrix giving where are the breakpoints each value of D (before choosing D): 
%       for every D, the vector to consider is rupt_1ZS(D,1:D)
% mu_1ZS : matrix giving the values at t_1, ..., t_n of the estimator ZS(D) for every D, 
%       the vector to consider is mu_1ZS(D,:)

%%%

n=length(Data);

%%% Complete infos with default values
% delta
if ~isfield(infos,'delta')
    infos.delta=1;
end

%%% Computes mh_ZS(D) for every D, and associated quantities

matD_ZhangSieg=cout_ZhangSieg( Data, infos.delta);
[contrast_ZhangSieg, rupt_1ZS]=prgdyn(matD_ZhangSieg,Dimmax);% nouvelles ruptures

compt=Inf*ones(1,Dimmax);
mu_1ZS=Inf*ones(Dimmax,n);
for d=1:Dimmax
    breaks=[0  rupt_1ZS( d, 1:d)];
    compte=diff(breaks);
    compt(d)=sum(log(compte));
    mu_1ZS( d, 1:compte(1))=sum(Data((breaks(1)+1):breaks(2)))/compte(1);
    for i=2:d
        mu_1ZS( d, (breaks(i)+1):breaks(i+1))=sum(Data((breaks(i)+1):breaks(i+1)))/compte(i)*ones(1,compte(i));
    end;
end;

%%% Selects the dimension D

crit2ZS = -(n-(1:Dimmax)+2)/2.*log(1+contrast_ZhangSieg)'+((1:Dimmax)-1)/2.*log(n*var(Data))-.5*compt+...
   log(n)*(1.5-(1:Dimmax))+gammaln((n-(1:Dimmax)+2)/2)-gammaln((n+1)/2);
D_ZS=find(crit2ZS==max(crit2ZS));

%%%

mu_ZS=mu_1ZS(D_ZS,:);
rupt_ZS=rupt_1ZS(D_ZS,1:D_ZS);


%%% Subfunctions used

function matD=cout_ZhangSieg( x, p)

n = length(x);
matD=Inf*ones(n,n);
x=reshape(x,1,n);

for i=1:n-p
    for j=i+p:n
        moy=mean(x(i:j));
        matD(i,j)=-(j-i+1)*(mean(x)-moy).^2/(n*var(x));
    end
end
