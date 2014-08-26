function matD=cout_Risque_emp( x, p)
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

n = length(x);
matD=Inf*ones(n,n);
x=reshape(x,1,n);

for i=1:n-p
    for j=i+p:n
        moy=mean(x(i:j));
        matD(i,j)=sum((x(i:j)-moy).^2)/n;
    end
end

