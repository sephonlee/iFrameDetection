function [Kd] = Kdiff(X1,X2,d,mode)

X1t1 = X1(:,1:size(X1,2)/2);
X1t2 = X1(:,size(X1,2)/2+1:end);
X2t1 = X2(:,1:size(X2,2)/2);
X2t2 = X2(:,size(X2,2)/2+1:end);

K11 = (kernelmatrix(mode{1},X1t1',X2t1',d(1)));
K22 = (kernelmatrix(mode{2},X1t2',X2t2',d(2)));
K12 = (kernelmatrix(mode{3},X1t1',X2t2',d(3)));

if length(d) == 4
    K21 = (kernelmatrix(mode{4},X1t2',X2t1',d(4)));
else
    K21 = K12'; % symmetric, same parameter
end

Kd = K11 + K22 - K12 - K21;
%       Kd = Kd+0.5.*eye(size(Kd,1));

end