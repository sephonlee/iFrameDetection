function [ D_hat, mu_hat, rupt_hat, rupt_D, s_D, Risk]=chpt1xx2VF( Data, A, p, nb_VF, D, infos)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input
% Data : vecteur de donnees
% A : methode de l'etape 1
%       A=1 : ERM
%       A=3 : Leave-p-out
% p : parametre caracterisant la methode 1 (V de V-fold / p de Lpo / parametre de la penalite Rademacher)
% nb_VF : nb de blocs de la VF (a l'etape 2)
% D : dimension maximale envisagee
% infos: structure contenant des informations annexes (mais utiles)
%       infos.delta=taille minimale d'un bloc - 1
%           delta=0 : on garde tous les modeles
%           delta=1 : on supprime les modeles avec un segment de taille 1
%           defaut: delta=1
%       [ infos.H_1_v infos.H_0_v infos.H_v] = sortie precalculee de som(n,p,1:n)
%           defaut : on la recalcule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Output
% D_hat : dimension finale selectionnee
% mu_hat : valeurs de l'estimateur final en les t_i
% rupt_hat : positions des ruptures de l'estimateur final
%       (avec la convention de eval_histo, les ruptures sont en:
%           t(1+rupt_hat(1:end-1))  
%           i.e., juste avant ces points)
% rupt_D : matrice contenant la meme information que rupt_hat (position des
%       ruptures de mh_1(dim)) pour toutes les valeurs de dim 
%       (avec des 0 pour completer)
% s_D : matrice contenant les valeurs de l'estimateur sh_(mh_1(dim)) 
%       en les t_i pour toutes les dimensions dim=1..D
% Risk : vecteur des valeurs de crit_2(dim) pour dim=1..D 
%       = estimations du risque de mh_1(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parametres
n=length(Data);
t=(1:n);
% Reg=zeros(size(Data));

% Complete infos avec les valeurs par defaut
% delta
if isfield(infos,'delta')
    delta=infos.delta;
else
    delta=1;
end
% si Lpo (A==3) , H_1_v H_0_v H_v
if ((A==3) && (~( isfield(infos,'H_1_v') && isfield(infos,'H_0_v') && isfield(infos,'H_v'))))
    [ infos.H_1_v infos.H_0_v infos.H_v] = som(n,p,1:n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 2: V-fold

% res=randperm(n);
aux=VFblocs(n,nb_VF);
[val,res]=sort(aux);
Z=Data(res);
indice=t(res);
risk=Inf*ones(nb_VF,D);

for v = 1:nb_VF
    Y_train =Z;
    Y_test = Z(find(val==v));
    Y_train(find(val==v))=[];
    n_train=length(Y_train);
    ind_train=indice;
    ind_train(find(val==v))=[];
    ind_test=indice(find(val==v));
    [t_train, perm_train]=sort(ind_train);
    Y_train=Y_train(perm_train);
    [t_test, perm_test]=sort(ind_test);
    Y_test=Y_test(perm_test);
%     Reg_train=Reg;
%     Reg_train(find(val==v))=[];
%     Reg_train=Reg_train(perm_train);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Etape 1
    contrast=Inf*ones(n_train,n_train);
    if (A==1)
        % Empirical contrast
        contrast = cout_Risque_emp( Y_train, delta);
        [contrast_min , rupt] = prgdyn( contrast, D);
    end;
    if (A==3)
        % LPO-risk
        contrast = cout_LPO( Y_train, p, delta, infos.H_1_v , infos.H_0_v , infos.H_v);
        [contrast_min , rupt] = prgdyn( contrast, D);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % suite de l'etape 2 
    for d = 1:D
        estim=Inf*ones(1,d);
        crit=Inf*ones(1,d);
        
        if d==1
            estim(1)=mean(Y_train);
            crit(1)=(sum(t_test<=t_train(rupt(1,1)))>=1)*sum((Y_test((t_test<=t_train(rupt(1,1))))-estim(1)).^2);
        else
            estim(1)=mean(Y_train(1:rupt(d,1)));
            crit(1)=(sum(t_test<=t_train(rupt(d,1)))>=1)*sum((Y_test((t_test<=t_train(rupt(d,1))))-estim(1)).^2);
            if d>2
                for k=2:d-1
                    estim(k)=mean(Y_train((rupt(d,k-1)+1):rupt(d,k)));
                    crit(k)=(sum(((t_test<=t_train(rupt(d,k))).*(t_test>t_train(rupt(d,k-1))))==1)>=1)*...
                        sum((Y_test(((t_test<=t_train(rupt(d,k))).*(t_test>t_train(rupt(d,k-1))))==1)-estim(k)).^2);
                end;
            end;
            estim(d)=mean(Y_train((rupt(d,d-1)+1):rupt(d,d)));
            crit(d)=(sum(t_test>t_train(rupt(d,d-1)))>=1)*sum((Y_test((t_test>t_train(rupt(d,d-1))))-estim(d)).^2);
        end;
        risk(v,d)=sum(crit(1:d))/length(t_test);
    end;
end;

Risk=mean(risk); % critere minimise pour choisir D = moyenne sur les V blocs pour chaque dimension d
[val,D_hat]=min(Risk);% D_hat : dimension selectionnee = celle qui minimise Risk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de l'estimateur bati sur toutes les donnees avec Etape 1
contrast2=Inf*ones(n,n);
if (A==1)
    contrast2 = cout_Risque_emp( Data, delta);
end;
if (A==3)
    contrast2 = cout_LPO( Data, p, delta, infos.H_1_v , infos.H_0_v , infos.H_v);
end;

[R_hat_D , rupt_D] = prgdyn( contrast2, D);
%%%%%%%%%%%%
% Calcul de l'erreur d'estimation pour chaque dimension
%%%%%%%%%%%%
% ell_D=Inf*ones(1,D);% vecteur des erreurs d'estimation pour chaque dimension
s_D=Inf*ones(D,n);% matrices des estimateurs pour chaque dimension
for d=1:D
    breaks=[0 rupt_D( d, 1:d)];
    compte=diff(breaks);
    s_D( d, 1:compte(1))=sum(Data((breaks(1)+1):breaks(2)))/compte(1);
    for i=2:d
        s_D( d, (breaks(i)+1):breaks(i+1))=sum(Data((breaks(i)+1):breaks(i+1)))/compte(i)*ones(1,compte(i));
    end;
end;
% ell_D=(mean((s_D-repmat(Reg,D,1)).^2,2))';
%%%%%%%%%%%%
% Calcul estimateur final
%%%%%%%%%%%%
% [ Risk_hat, blabla]=min(R_hat_D);
R_hat_D=R_hat_D';
% Risk_hat=min(R_hat_D);
mu_hat=s_D(D_hat,1:n);
% ell_hat=ell_D(D_hat);
rupt_hat=rupt_D( D_hat, 1:D_hat);


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


function contrast = cout_LPO( Data, p, delta, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce programme permet le calcul de la matrice de cout a implementer.
%
% Data: est le vecteur de donnees
% p: est le p du LPO
% delta
% argument optionnel: precalcul du resultat de som(n,p,1:n) i.e., 3
%   vecteurs dans l'ordre H_1_v, H_0_v et H_v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
n=length(Data);

if size(varargin,2)>2
    H_1_v=varargin{1};
    H_0_v=varargin{2};
    H_v=varargin{3};
else
    [ H_1_v, H_0_v, H_v]=som(n,p,1:n);
end

contrast=Inf*ones(n,n);
for i=1:n
    for j=i:n
        % LPO contraste
%         [ H_1, H_0, H]=som(n,p,j-i+1);
        compt=j-i+1;

        H_1=H_1_v(compt);
        H_0=H_0_v(compt);
        H=H_v(compt);
        
        
        if compt<=max(1,delta)
            contrast(i,j)=Inf;
        else
            if compt<p
               q0 = exp(gammaln(p+1)-gammaln(p-compt+1)-(gammaln(n+1)-gammaln(n-compt+1)));

               cste = 1-exp(gammaln(n-compt+1)-gammaln(n-p+1)-gammaln(p-compt+1)-(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1)));
                contrast(i,j)=sum(Data(i:j).^2)*(1/n-1/p*q0)+...
                    1/p*sum(Data(i:j).^2)*(H_1-H_0/compt)-...
                    2/p*((sum(Data(i:j)))^2-sum(Data(i:j).^2))*(H_0/(compt-1)-H/compt/(compt-1));
                if (compt>2)
                    contrast(i,j)=contrast(i,j)+...
                        ((compt-2)*(sum(Data(i:j)))^2-(compt-2)*sum(Data(i:j).^2))/p*((compt+1)*H_0-H-compt*H_1)/compt/(compt-1)/(compt-2);
                end;
                contrast(i,j)=contrast(i,j)/cste;
            else
                contrast(i,j)=sum(Data(i:j).^2)/n+...
                    1/p*sum(Data(i:j).^2)*(H_1-H_0/compt)-...
                    2/p*((sum(Data(i:j)))^2-sum(Data(i:j).^2))*(H_0/(compt-1)-H/compt/(compt-1));
                if (compt>2)
                    contrast(i,j)=contrast(i,j)+...
                        ((compt-2)*(sum(Data(i:j)))^2-(compt-2)*sum(Data(i:j).^2))/p*((compt+1)*H_0-H-compt*H_1)/compt/(compt-1)/(compt-2);
                end;
            end;
        end;
             
    end;
end;


function [ H_1, H_0, H]=som(n,p,compt)

D=length(compt);
H_1=Inf*ones(1,D);
H_0=Inf*ones(1,D);
H=Inf*ones(1,D);
for ind=1:D
    q=zeros(1,compt(ind));

    r=max(1,compt(ind)-p):min(compt(ind),n-p);

    q(r)=exp(gammaln(n-p+1)-gammaln(r+1)-gammaln(n-p-r+1)-...
                (gammaln(n+1)-gammaln(compt(ind)+1)-gammaln(n-compt(ind)+1))+...
                gammaln(p+1)-gammaln(compt(ind)-r+1)-gammaln(p-compt(ind)+r+1));
 
    inv=1./r;
    H_1(ind)=sum(inv.*q(r));
    H_0(ind)=sum(q);
    H(ind)=sum(r.*q(r));
end;



