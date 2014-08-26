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