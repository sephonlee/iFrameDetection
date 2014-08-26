% Test Accuracy
% clear
% close all
% root = '/media/MITCH/Data/UNSU_CD_KKM/KernelDiff/public';
% ID = ['Gr_isl_' num2str(OO.NtrnPix) '/'];
% DIR = [root '/RESULTS/SETS/'];
% addpath(genpath(root))
% cd(root)
%%%%% GLOBAL PARAMETERS
% OO.Nclust = 2;
% OO.trnSetsSizes = [100];
% OO.plots = 'off'; % 'off'
% OO.Nsets = [1 2]; % # experiments, can be a vector with specific indexes, e.g. [2 4 5] ocio a read exp
% OO.kernels = {'RBF'};%,};%{'lin'},
% OO.methods = {'diff'};%
% OO.mu = 1;

%%%% LOAD LABELS to test            
load('./imm/Greece_Island_ROI.mat');
labels_change = ROI;

for jj = 1:length(OO.methods)
    OO.method = OO.methods{jj};
    for kk = 1:length(OO.kernels)
        OO.ktype = OO.kernels{kk};
        for ii = OO.trnSetsSizes
            OO.NtrnPix = ii; % == the training set to use
            if strcmp(OO.method,'kernelDiff')
                OO.ktype = {OO.kernels{kk},OO.kernels{kk},OO.kernels{kk},OO.kernels{kk}};
                %               clso   OO.ktype = {{'lin'},{'lin'},{'lin'},{'lin'}};
                OO.ktypedisp = [OO.ktype{1} OO.ktype{3}];
            else
                %                 OO.ktype = 'lin';
                OO.ktypedisp = OO.ktype;
            end
            
            %%%%% ID SUBSET (Subfolder in which it saves results)
            
            load([DIR ID 'EXP_' OO.method '_ktype_' OO.ktypedisp '_' num2str(OO.NtrnPix) '_Map_sum.mat'],'Map_sum');

            
            
            % T_map = data2column(labels_C_NC);
            T_map = data2column(labels_change);
            T_lab = T_map(T_map ~= 0);
            P_maps = data2column(Map_sum);
            P_lab = P_maps(T_map ~= 0,:);
            Assessment.Kappa = zeros(1,length(OO.Nsets));
            Assessment.OA = zeros(1,length(OO.Nsets));
            for i = OO.Nsets
                load([DIR ID 'EXP_' OO.method '_ktype_' OO.ktypedisp ...
                    '_' num2str(OO.NtrnPix) 'trn' num2str(i) '.mat'],'dist');
                [TEST] = assessment(T_lab,P_lab(:,i)+1,'class');
                Assessment.Kappa(i) = TEST.Kappa;
                Assessment.OA(i) = TEST.OA;
                
                if  mean(dist(T_map==2,1)-dist(T_map==2,2)) > mean(dist(T_map==1,1)-dist(T_map==1,2))
                    dista = dist(T_map~=0,1)-dist(T_map~=0,2);
                else dista = dist(T_map~=0,2)-dist(T_map~=0,1);
                end
                
                %                 [AUC,tpr,fpr]=roccurve(P_lab(:,i),T_lab-1,'noPlot');
                [tpr,fpr,t,AUC]=perfcurve(T_lab,dista,2);
                %                 figure; plot(tpr,fpr,'Color','black'), hold on
                %                 plot(tpr2,fpr2,'Color','red')
                Assessment.AUC(i) = [AUC];
                Assessment.tpr{i} = tpr;
                Assessment.fpr{i} = fpr;
                Assessment.AdjustedRand(i) = adjrand(P_lab(:,i)+1,T_lab);
                clear TEST
                fprintf('set%i',i)
            end
            Assessment.mAdjustedRand = mean(Assessment.AdjustedRand);
            Assessment.mKappa = mean(Assessment.Kappa);
            Assessment.mOA = mean(Assessment.OA);
            Assessment.mAUC = mean(Assessment.AUC);
            Assessment.sAdjustedRand = std(Assessment.AdjustedRand);
            Assessment.sKappa = std(Assessment.Kappa);
            Assessment.sOA = std(Assessment.OA);
            Assessment.sAUC = std(Assessment.AUC);
            
            fprintf('.kernel_%s.method_%s\n',OO.ktypedisp,OO.method)
            
            save([DIR ID 'EXP_' OO.method '_ktype_' OO.ktypedisp '_' num2str(OO.NtrnPix) '_Map_sum.mat'],'Assessment','-append');
        end
        
        
    end
end
