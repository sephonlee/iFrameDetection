function [WholeSet VAL fitGMM] = initialization_KKM(imm1,imm2,OO,savedir)
Nsets = OO.Nsets;
trnSetSizes = OO.trnSetsSizes;

if size(imm1,3)~= 1
    imm1 = data2column(imm1);
    imm2 = data2column(imm2);
end

if ~strcmp(class(imm2),'double')
    imm1 = double(imm1);
    imm2 = double(imm2);
end

WholeSet = [imm1 imm2];

diff = imm2-imm1;
VAL= [sqrt(sum(diff.^2,2))]; 
% the first time you run automatic intialization, check the magnitude!
% hist(VAL,50)
% figure; imagesc(reshape(VAL,OO.OriginalSize(1),OO.OriginalSize(2))); axis image
fitGMM = gmdistribution.fit(VAL,2,'Start','randSample','CovType','diagonal');%,'Options',options
mu1 = min(fitGMM.mu);
mu2 = max(fitGMM.mu);
sig1 = min(fitGMM.Sigma);
sig2 = max(fitGMM.Sigma);

% CHANGE HERE IF EVERYTHING'S WRONG: This initialization is sometimes difficult 

thresh = [0 mu1+0.5 mu2-0.5 max(VAL)];

ch = find(VAL(:,1)>=thresh(3) & VAL(:,1)<=thresh(4));
n_ch = find(VAL(:,1)>thresh(1) & VAL(:,1)<thresh(2));

for j = 1:length(trnSetSizes)
    no_trnpts = trnSetSizes(j);
    for i = Nsets
        rans_c = randperm(length(ch))';
        rans_nch = randperm(length(n_ch))';
        cord = [ch(rans_c(1:no_trnpts/2,:)); n_ch(rans_nch(1:no_trnpts/2,:))];
        trainSet = [WholeSet(ch(rans_c(1:no_trnpts/2)),:); WholeSet(n_ch(rans_nch(1:no_trnpts/2)),:)];
        
        mean_c = median(trainSet(1:no_trnpts/2,:));
        mean_nc = median(trainSet(no_trnpts/2+1:end,:));
        
        [temp,I_c] = min(sum((trainSet(1:no_trnpts/2,:) - repmat(mean_c,no_trnpts/2,1))));
        [temp,I_nc] = min(sum((trainSet(no_trnpts/2+1:end,:) - repmat(mean_nc,no_trnpts/2,1))));
        I_nc = I_nc + no_trnpts/2;
        
        centers = [WholeSet(I_c,:);WholeSet(I_nc,:)];
        
        if ~isdir(savedir)
            mkdir(savedir)
        end
        save([savedir '/trainSet' num2str(i) '.mat'],'trainSet','centers','cord')
    end
end


