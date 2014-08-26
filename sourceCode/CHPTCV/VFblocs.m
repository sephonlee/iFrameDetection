function vect=VFblocs(n,V)
% fonction retournant un decoupage en V blocs de 1:n
% de maniere aussi reguliere que possible
% sous la forme d'un vecteur de taille n dont les elements sont
% des entiers compris entre 1 et V (= le num du bloc ou ils sont)
% pour retrouver la liste des indices du bloc j, faire find(vect==j)
vect=zeros(1,n);
t=floor(n/V);
r=n-V*t;
if r>0
    s=floor((r+1)*rand(1));
    l=randperm(V);
	if s>0
        vect(1:s)=l(1:s);
	end%s>0
    if (r-s)>0
        vect((n-r+s+1):n)=l((s+1):r);
	end%(r-s)>0
else
    s=0;
end%r>0
%
for j=0:(t-1)
    vect((s+1+j*V):(s+(j+1)*V))=randperm(V);
end%j=1:V