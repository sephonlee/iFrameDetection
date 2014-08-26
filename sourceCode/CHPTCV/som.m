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