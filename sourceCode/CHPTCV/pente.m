% *************************************************
% ***    calcul du alpha "optimal"      	***
% *************************************************
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
function [k_max,k_seuil,alpha_max,alpha_seuil]=pente(J,pen,seuil)
K=length(J);
% pen=(1:K)'.*(c2+c1*log(n./(1:K)'));
k=1;
kv=[];
dv=[];
pv=[];
dmax=1;
while k<K
	pk=(J(k+1:K)-J(k))./(pen(k)-pen(k+1:K));
	[pm,dm]=max(pk);
	kv=[kv k];
	dv=[dv dm];
	pv=[pv pm];
	if dm>dmax;  dmax=dm; k_max=k; pmax=pm;  end;
	k=k+dm;
end;

pv=[pv 0];
kv=[kv K];
dv=diff(kv);
[dmax,rmax]=max(dv);
rt=find(dv==dmax);
pmax=pv(rt(end));
alpha_max=2*pmax;
rt=find(alpha_max>=pv);
k_max=kv(rt(1));

%%% Ajout de Sylvain: Kmin determinee par le moment ou la dimension passe
%%% sous seuil
p_seuil=min(pv(kv<=seuil));
alpha_seuil=2*p_seuil;
rt=find(alpha_seuil>=pv);
k_seuil=kv(rt(1));
