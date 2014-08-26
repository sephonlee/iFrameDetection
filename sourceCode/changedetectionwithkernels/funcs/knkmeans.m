function [label, cent_ind, DPC, DBC] = knkmeans(K,m)

% K: kernel matrix
% m: k (integer) or label (1 x n, 1<=label(i)<=k)

% reference: [1] Kernel Methods for Pattern Analysis
% by John Shawe-Taylor, Nello Cristianini

% Written by Mo Chen (mochen@ie.cuhk.edu.hk). March 2009. 
% Modified by Michele Volpi 2011, For: Unsuopervised change detection with
% kernels, IEEE GRSL, 2012.

n = size(K,1);
if max(size(m)) == 1
    k = m;
%     label = randi(k,1,n);
    label = ceil(k*rand(1,n));
elseif size(m,1) == 1 && size(m,2) == n
    k = max(m);
    label = m;
else
    error('ERROR: m is not valid.');
end
conv = 1; % distances_new - distances_old
conv_thresh = 1000; % convergence threshold
last = 0;
S = repmat((1:k)',1,n);

while any(label ~= last) && conv < conv_thresh
    [temp1,temp2,label] = unique(label); % remove empty clusters
    L = double(bsxfun(@eq,S,label)); % 1-0 cluster indicator 
    E = bsxfun(@rdivide,L,sum(L,2)); % 1/n_clust 
    T = E*K; % K / 1/n_clust (T(c,:) = dist x_i - cluster center)
    Z = repmat(diag(K),1,k)' + repmat(diag(T*E'),1,n)-2*T;  % dist(tutti in cluster c) -2*(dist dal centro per cluster c)
    last = label; % check update rule
    [temp3, label] = min(Z); % min fra i dist al centro (dist dal cluster = linea)
    conv = conv+1;
end

[temp cent_ind] = min(Z,[],2);

DPC = zeros(1,m);
for j = 1:m
idx = label == j;
DPC(j) = sum(Z(j,idx),2)'; % trace(K(idx,idx))+
end

% D_fra i clust in H:
if m == 2 %only for 2 clusters, generalize
    id1 = label == 1;
    id2 = label == 2;
    ell1 = sum(id1);
    ell2 = sum(id2);
a = sum(sum(K(id1,id1)))/ell1^2;
b = sum(sum(K(id2,id2)))/ell2^2;
c = sum(sum(K(id2,id1)))/(ell1*ell2);
DBC = a+b-2*c;

else
    DBC = 0;
end




