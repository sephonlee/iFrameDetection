function [D_PML, mu_PML, rupt_PML, crit2PML, rupt_1PML, mu_1PML, D_PMLmax, mu_PMLmax, rupt_PMLmax, crit2PMLmax]=proc_PML(Data, Dimmax, infos)
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
% Procedure [PML]
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
% D_PML : dimension selected by the procedure [PML] ('PML' in the paper)
% mu_PML : values of the final estimator at t_1, ..., t_n
% rupt_PML : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2PML : vector of estimated values of the risk (with the PML criterion) 
%           for all values of the dimension between 1 and Dimmax 
%
% rupt_1PML : matrix giving where are the breakpoints each value of D (before choosing D): 
%       for every D, the vector to consider is rupt_1PML(D,1:D)
% mu_1PML : matrix giving the values at t_1, ..., t_n of the estimator PML(D) for every D, 
%       the vector to consider is mu_1PML(D,:)
%
% D_PMLmax : dimension selected by the procedure [PML] with a different version of the slope heuristics ('max jump' instead of 'threshold')
% mu_PMLmax : values of the final estimator at t_1, ..., t_n
% rupt_PMLmax : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2PMLmax : vector of estimated values of the risk (with the PML criterion) 
%           for all values of the dimension between 1 and Dimmax 
%

%%%

n=length(Data);

%%% Complete infos with default values
% delta
if ~isfield(infos,'delta')
    infos.delta=1;
end

%%% Computes mh_PML(D) for every D, and associated quantities

matD_Pml=cout_Pml( Data, infos.delta);% delta=1 car base sur un estimateur de la variance dans chaque case
[contrast_Pml, rupt_1PML]=prgdyn(matD_Pml,Dimmax);

mu_1PML=Inf*ones(Dimmax,n);
for d=1:Dimmax
    breaks=[0  rupt_1PML( d, 1:d)];
    compte=diff(breaks);
    mu_1PML( d, 1:compte(1))=sum(Data((breaks(1)+1):breaks(2)))/compte(1);
    for i=2:d
        mu_1PML( d, (breaks(i)+1):breaks(i+1))=sum(Data((breaks(i)+1):breaks(i+1)))/compte(i)*ones(1,compte(i));
    end;
end;

%%% Selects the dimension D

pen_lin=(1:Dimmax);
[D_PMLmax,D_PML,alpha_max_Pml,alpha_seuil_Pml]=pente(contrast_Pml',pen_lin,floor(n/log(n)));
crit2PMLmax=contrast_Pml'+alpha_max_Pml*pen_lin;
crit2PML=contrast_Pml'+alpha_seuil_Pml*pen_lin;

%%%

mu_PML=mu_1PML(D_PML,:);
rupt_PML=rupt_1PML(D_PML,1:D_PML);

mu_PMLmax=mu_1PML(D_PMLmax,:);
rupt_PMLmax=rupt_1PML(D_PMLmax,1:D_PMLmax);


%%% Subfunctions used

function matD=cout_Pml( x, p)

n = length(x);
matD=Inf*ones(n,n);
x=reshape(x,1,n);

for i=1:n-p
    for j=i+p:n
        moy=mean(x(i:j));
        matD(i,j)=(j-i+1)*log(mean((x(i:j)-moy).^2));
    end
end
