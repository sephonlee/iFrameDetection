function [D_ERMBMest, mu_ERMBMest, rupt_ERMBMest, crit2BMest_1ERM, D_ERMBMthr, mu_ERMBMthr, rupt_ERMBMthr, crit2BMthr_1ERM, D_ERMBMmax, mu_ERMBMmax, rupt_ERMBMmax, crit2BMmax_1ERM, rupt_1ERM, mu_1ERM]=proc_ERMBM(Data, Dimmax, infos)
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
% Procedure [ERM,BM]
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
%       infos.threshold= threshold used in the slope heuristics algorithm
%           (see Section 4 of Supplementary material for details)
%           defaut: threshold=min(floor(n/log(n)),floor(0.75*Dimmax))
%
% Output:
%
% D_ERMBMest : dimension selected by the procedure [ERM,BM] with sigma estimated ('BM' in the paper)
% mu_ERMBMest : values of the final estimator at t_1, ..., t_n
% rupt_ERMBMest : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2BMest_1ERM : vector of estimated values of the risk (with the Birge-Massart criterion) 
%           for all values of the dimension between 1 and Dimmax 
%
% D_ERMBMthr : dimension selected by the procedure [ERM,BM] with the slope
%       heuristics + threshold ('BM_{thresh.}' in Section 4 of the Supplementary material)
% mu_ERMBMthr : values of the final estimator at t_1, ..., t_n
% rupt_ERMBMthr : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2BMthr_1ERM : vector of estimated values of the risk (with the Birge-Massart criterion) 
%           for all values of the dimension between 1 and Dimmax 
%
% D_ERMBMmax : dimension selected by the procedure [ERM,BM] with the slope
%       heuristics + maximal jump ('BM_{max.jump}' in Section 4 of the Supplementary material)
% mu_ERMBMmax : values of the final estimator at t_1, ..., t_n
% rupt_ERMBMmax : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2BMmax_1ERM : vector of estimated values of the risk (with the Birge-Massart criterion) 
%           for all values of the dimension between 1 and Dimmax 
%
% rupt_1ERM : matrix giving where are the breakpoints each value of D (before using BM for choosing D): 
%       for every D, the vector to consider is rupt_1ERM(D,1:D)
% mu_1ERM : matrix giving the values at t_1, ..., t_n of the estimator ERM(D) for every D, 
%       the vector to consider is mu_1ERM(D,:)

%%%

n=length(Data);

%%% Complete infos with default values
% delta
if ~isfield(infos,'delta')
    infos.delta=1;
end
% threshold
if ~isfield(infos,'threshold')
    infos.threshold=min(floor(n/log(n)),floor(0.75*Dimmax));
end

%%% Computes mh_ERM(D) for every D, and associated quantities

matD=cout_Risque_emp( Data, infos.delta);
[contrast,rupt_1ERM]=prgdyn(matD,Dimmax);

mu_1ERM=Inf*ones(Dimmax,n);% matrices des estimateurs pour chaque dimension
for d=1:Dimmax
    breaks=[0 rupt_1ERM( d, 1:d)];
    compte=diff(breaks);
    mu_1ERM( d, 1:compte(1))=sum(Data((breaks(1)+1):breaks(2)))/compte(1);
    for i=2:d
        mu_1ERM( d, (breaks(i)+1):breaks(i+1))=sum(Data((breaks(i)+1):breaks(i+1)))/compte(i)*ones(1,compte(i));
    end;
end;

%%% Selects the dimension D (several alternative methods possible here)

%penalite lineaire avec le log(n/D) + heuristique de pente
pen_lin_log=(1:Dimmax).*(5+2*log(n./(1:Dimmax)));
[D_ERMBMmax,D_ERMBMthr,alpha_max_lin_log,alpha_seuil_lin_log]=pente(contrast',pen_lin_log,infos.threshold);
crit2BMmax_1ERM=contrast'+alpha_max_lin_log*pen_lin_log;
crit2BMthr_1ERM=contrast'+alpha_seuil_lin_log*pen_lin_log;

%penalite lineaire avec le log(n/D) + sigma inconnu (estimateur de Baraud 2000)
sigma_estime=sqrt(estim_variance((1:n),Data));
alpha_sig_estime=(sigma_estime^2)/n;
crit2BMest_1ERM=contrast'+alpha_sig_estime*pen_lin_log;
D_ERMBMest=find(crit2BMest_1ERM==min(crit2BMest_1ERM),1,'first');
%


%%%

mu_ERMBMest=mu_1ERM(D_ERMBMest,:);
rupt_ERMBMest=rupt_1ERM(D_ERMBMest,1:D_ERMBMest);

mu_ERMBMmax=mu_1ERM(D_ERMBMmax,:);
rupt_ERMBMmax=rupt_1ERM(D_ERMBMmax,1:D_ERMBMmax);

mu_ERMBMthr=mu_1ERM(D_ERMBMthr,:);
rupt_ERMBMthr=rupt_1ERM(D_ERMBMthr,1:D_ERMBMthr);


%%% Subfunctions used

function v=estim_variance(X,Y)
% This function returns a classical estimate of the variance of an
% homoscedastic sequence, assuming implicitly that the mean of the function
% is smooth, i.e., varies slowly with X
n=numel(Y);
[b,ix]=sort(X);
c=Y(ix);
if n<=1
    echo('Warning: cannot estimate the variance with only one value');
    v=NaN;
else
    if mod(n,2)~=0
        c(n-1)=(c(n-1)+c(n))/2;
        c=c(1:(n-1));
        n=n-1;
    end
    v=sum((c(1:2:(n-1))-c(2:2:n)).^2)/n;
end

 
