function [D_ERMVF, mu_ERMVF, rupt_ERMVF, crit2VF_1ERM, rupt_1ERM, mu_1ERM]=proc_ERMVF(Data, nb_VF, Dimmax, infos)
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
% Procedure [ERM,VF]
%
% Input: 
% Data : vector of observations Y_1, ..., Y_n
% nb_VF : number of folds (V)
% Dimmax : maximal dimension considered
% infos: structure with additional informations
%       infos.delta= (minimal size of a segment)-1
%           delta=0 : all models are kept
%           delta=1 : models with at least one segment of the form
%           [t_i,t_{i+1}) are removed
%           defaut: delta=1
%
% Output:
% D_ERMVF : dimension selected by the procedure
% mu_ERMVF : values of the final estimator at t_1, ..., t_n
% rupt_ERMVF : values of the index i corresponding to breakpoints of the final estimator at t_i 
%           (with the convention used for s and sigma in Regsig_gener.m, plus a fictive breakpoint at t_{n+1})
% crit2VF_1ERM : vector of estimated values of the risk (via V-Fold cross-validation) for all values of the dimension between 1 and Dimmax 
% rupt_1ERM : matrix giving where are the breakpoints each value of D (before using VFCV for choosing D): 
%       for every D, the vector to consider is rupt_1ERM(D,1:D)
% mu_1ERM : matrix giving the values at t_1, ..., t_n of the estimator ERM(D) for every D, 
%       the vector to consider is mu_1ERM(D,:)

[ D_ERMVF, mu_ERMVF, rupt_ERMVF, rupt_1ERM, mu_1ERM, crit2VF_1ERM]=chpt1xx2VF( Data, 1, 0, nb_VF, Dimmax, infos);

 
