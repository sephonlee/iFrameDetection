% clear
% close all
% root = '/media/MITCH/Data/UNSU_CD_KKM/KernelDiff/';

root = '//psf/Home/Desktop/Research/PNNL/sourceCode/changedetectionwithkernels/funcs';
addpath(genpath(root))
cd(root)
%%%%%%%%%%% SINGLE KERNEL KERNEL K-MEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GLOBAL PARAMETERS
OO.Nclust = 2;
OO.trnSetsSizes = [100];%
OO.plots = 'off'; % 'off' | 'on'
OO.Nsets = [1:5];% # experiments, can be a vector with specific indexes, e.g. [2 4 5] ocio a read exp
OO.kernels = {{'RBF'}};% {'lin'},
OO.methods = {{'kernelDiff'}};%{'stack'},{'kernelDiff'},
% OO.methods = {{'diff'}};%{'diff'},{'stack'},
OO.mu = 1;
for jj = 1:length(OO.methods)
    OO.method = OO.methods{jj}{1};
    for kk = 1:length(OO.kernels)
        OO.ktype = OO.kernels{kk}{1};
        for ii = OO.trnSetsSizes
            OO.NtrnPix = ii; % == the training set to use
            if strcmp(OO.method,'kernelDiff')
                OO.ktype = {{OO.kernels{kk}{1}},{OO.kernels{kk}{1}},{OO.kernels{kk}{1}},{OO.kernels{kk}{1}}};
                % OO.ktype = {{'lin'},{'lin'},{'lin'},{'lin'}};
                OO.ktypedisp = [OO.ktype{1}{1} OO.ktype{3}{1}];
            else
                % OO.ktype = 'lin';
                OO.ktypedisp = OO.ktype{1}{1};
            end
            close all
            
            DIR = [root '/RESULTS/SETS/'];
            %%%%% LOAD IMAGE
            % load('F:/MITCH/DATA/UNSU_CD_KKM/KernelDiff/imm/VHR_QB_Bruttisellen_SMALL.mat')
            % load('./imm/DFC_subset.mat')
            %             load([ root '/imm/VHR_QB_Bruttisellen_SMALL.mat'])
            %             load([root '/imm/Greece_forestfires.mat'])
            load('./imm/Greece_island.mat');
            t1(:,:,6) = [];
            t2(:,:,6) = [];
            % load([root '/imm/Greece_forestfires.mat']);
            % load('F:/MITCH/DATA/UNSU_CD_KKM/KernelDiff/imm/set_vicino_zugo.mat')
            %%%%% ID SUBSET (Subfolder in which it saves results)
            % ID = ['Bruttisellen_long_NtrnPix_' num2str(OO.NtrnPix) '/'];
            %             ID = ['GreeceBur_NtrnPix_' num2str(OO.NtrnPix) '/'];
            ID = ['GreeceIsl_NtrnPix_' num2str(OO.NtrnPix) '/'];
            % ID = ['spotDFC_NtrnPix_' num2str(OO.NtrnPix) '/'];
            %             ID = ['TEST_Bruttisellen_small_NtrnPix_' num2str(OO.NtrnPix) '/'];
            %             ID = ['WOTIR_GreeceBur_NtrnPix_' num2str(OO.NtrnPix) '/'];
            % ID = ['BruttiMS_NtrnPix_' num2str(OO.NtrnPix) '/'];
            %             ID = ['TEST_VHRSMALL_' num2str(OO.NtrnPix) '/'];
            [o1 o2 o3] = size(t1);
            OO.OriginalSize = [o1 o2 o3];
            %%%%% Create the complete set to cluster
            d1 = data2column(standardimg(double(t1),'zscores'));
            d2 = data2column(standardimg(double(t2),'zscores'));
            
            switch OO.method
                case 'stack'
                    WholeSet = [d1 d2]; % FOR THE STACK APPROACH
                    [s1 s2] = size(WholeSet);
                    OO.ColSize = [s1 s2];
                case 'diff'
                    WholeSet = [d2-d1]; % FOR THE DIFFERENCE IMAGE CLUSTERING
                    [s1 s2] = size(WholeSet);
                    OO.ColSize = [s1 s2];
                case 'kernelDiff'
                    WholeSet = [d1 d2]; % FOR THE DIFFERENCE IN THE FEATURE SPACE CLUSTERING
                    [s1 s2] = size(WholeSet);
                    OO.ColSize = [s1 s2];
            end
            
            for i = OO.Nsets
                disp(i)
                fprintf(['Single Kernel (' OO.ktypedisp ') run #: %i\n'],i);
                st = sprintf('load([DIR ID ''trainSet%i.mat''],''trainSet'',''centers'',''cord'')',i);
                eval(st)
                centers(:,[6 13]) = [];
                trainSet(:,[6 13]) = [];
                % eval(sprintf('centers = centers%i;',i))
                % eval(sprintf('trainSet = trainSet%i;',i))
                % eval(sprintf('cord = cord%i;',i))
                % eval(sprintf('clear centers%i trainSet%i cord%i cord',i,i,i)) % ERASE CORD IF NEEDED
                %                 trainSet(:,[6 13]) = []; % se rinouverre termico da landsat
                % clear centers trainSet cord % ERASE CORD IF NEEDED
                %%%%% Prepare Centers / TrainSet According to the method:
                if strcmp(OO.method,'diff')
                    % centers in difference coordinates
                    centers = centers(:,OO.ColSize(2)+1:OO.ColSize(2)*2)-centers(:,1:OO.ColSize(2));
                    % centers in difference coordinates
                    trainSet = trainSet(:,OO.ColSize(2)+1:OO.ColSize(2)*2)-trainSet(:,1:OO.ColSize(2));
                end
                %%%%% Inital guess on the sigma parameter (IF NEEDED)
                if iscell(OO.ktype) && strcmp(OO.ktype{1}{1},'RBF')
                    fprintf('RBF sigma Guess: ');
                    randind = randperm(OO.ColSize(1))';
                    OO.kpar = 1/2*(median(median(L2_distance(WholeSet(randind(1:3000),1:OO.ColSize(2)/2),WholeSet(randind(1:3000),1:OO.ColSize(2)/2)))) ...
                        + median(median(L2_distance(WholeSet(randind(1:3000),OO.ColSize(2)/2+1:OO.ColSize(2)),WholeSet(randind(1:3000),OO.ColSize(2)/2+1:OO.ColSize(2))))));
                    % OO.kpar =
                    % median(median(L2_distance(trainSet(size(trainSet,1)/2+1:end,:),trainSet(size(trainSet,1)/2+1:end,:))));
                    % OO.kpar = median(median(L2_distance(trainSet,trainSet)));
                    % OO.kpar = median(median(L2_distance(trainSet(size(trainSet,1)/2+1:end,:),trainSet(size(trainSet,1)/2+1:end,:))));
                    fprintf([num2str(OO.kpar) '\n']);
                elseif ~iscell(OO.ktype) && (strcmp(OO.ktype,'RBF') | strcmp(OO.ktype,'sam'))
                    fprintf('RBF sigma Guess: ');
                    randind = randperm(OO.ColSize(1))';
                    OO.kpar = median(median(L2_distance(WholeSet(randind(1:3000),:),WholeSet(randind(1:3000),:))));
                    % OO.kpar = median(median(L2_distance(trainSet(size(trainSet,1)/2+1:end,:),trainSet(size(trainSet,1)/2+1:end,:))));
                    % OO.kpar = median(median(L2_distance(trainSet,trainSet)));
                    fprintf([num2str(OO.kpar) '\n']);
                elseif strcmp(OO.ktype,'poly')
                    OO.kpar = [1 2 3];
                else
                    OO.kpar = 1;
                end
                clear t1 t2 d1 d2 randind NtrnPix o1 o2 o3 s1 s2% clean up
                %%%%% Make the search among the candidate parameter (par vector -> pars)
                switch OO.method
                    case {'diff','stack'}
                        if ~strcmp(OO.ktype,'lin')
                            % pars = 0.5*OO.kpar:0.2:5*OO.kpar;
                            pars = 0.1:0.1:10;
                        else
                            pars = 1;
                        end
                        %%% initialize tab_errors
                        tab_errors = zeros(length(pars),8);
                        lab_mem = [ones(size(trainSet,1)/2,1); 2*ones(size(trainSet,1)/2,1)]';
                        for p1 = pars
                            fprintf(['PARAMETER SEARCH: ' num2str(min(pars)) ' to ' num2str(max(pars)) ...
                                ' by ' num2str(0.1) '....now testing: ' num2str(p1) '\n']);
                            %%%%% BUILD KERNEL
                            K = kernelmatrix(OO.ktype,trainSet',trainSet',p1);
                            % K = normalise(K);
                            % K = (p1)*kernelmatrix('RBF',trainSet(:,1:4)',trainSet(:,1:4)',prior1) + ...
                            %     (1-p1)*kernelmatrix('RBF',trainSet(:,5:8)',trainSet(:,5:8)',prior2);
                            %%%%% Train The Kernel k-Means
                            [lab_trn, trn_centers, DPC, DBC] = knkmeans(K,OO.Nclust);
                            id1 = (lab_trn == 1)';
                            id2 = (lab_trn == 2)';
                            ell1 = sum(id1);
                            ell2 = sum(id2);
                            % DPCs = DPC(1,1)./ell1 + DPC(1,2)./ell2; % Normalize DPC
                            DPCs = (1/(OO.NtrnPix))*(sum(abs(DPC(:)))); % Normalize DPC
                            % cost =((1-OO.mu)*DPCs-OO.mu*DBC);   %
                            % cost = DPCs/DBC;   %
                            cost = (DPCs-DBC);   %
                            %                             cost = DPCs/DBC^2;
                            tab_errors((pars == p1)',:) = [p1 trn_centers' ell1 ell2 DPCs DBC^2 cost];
                            clear id1 id2 ell1 ell2 DPCs DBC cost K
                        end
                    case 'kernelDiff'
                        % pars1 = OO.kpar;
                        % pars2 = 1;
                        if strcmp(OO.ktypedisp,'linlin')
                            pars1 = 1;
                            pars2 = 1;
                        else
                            %                                                         pars1 = OO.kpar;
                            %                                                         pars2 = 2*OO.kpar;
                            pars1 = 1:0.1:3; %
                            pars2 = 5:0.1:6; %
                        end
                        %%% initialize tab_errors
                        tab_errors = zeros(length(pars1)*length(pars2),11);
                        count = 1;
                        for p1 = pars1
                            for p2 = pars2
                                pp = [p1 p1 p2 p2];
                                fprintf(['PARAMETER SEARCH: now testing: ' num2str(pp(1)) ' ' num2str(pp(2)) ...
                                    ' ' num2str(pp(3)) ' ' num2str(pp(4)) '\n']);
                                %%%%% BUILD KERNEL
                                K = Kdiff(trainSet,trainSet,pp,OO.ktype);
                                %                                 K = centering(K);
                                % K = normalise(K);
                                % ps = [1*ones(OO.trnSetsSizes/2,1); 2*ones(OO.trnSetsSizes/2,1)]';
                                [lab_trn, trn_centers, DPC, DBC] = knkmeans(K,OO.Nclust);
                                id1 = (lab_trn == 1)';
                                id2 = (lab_trn == 2)';
                                ell1 = sum(id1);
                                ell2 = sum(id2);
                                % DPCs = DPC(1,1)./ell1 + DPC(1,2)./ell2; % Normalize DPC
                                % DPCs = abs(1/OO.NtrnPix.*(sum(DPC(:)))); % Normalize DPC
                                DPCs = 1/OO.NtrnPix*(sum(abs(DPC(:)))); % Normalize DPC
                                % cost = ((1-OO.mu)*DPCs-OO.mu*DBC);   %
                                % cost = DPCs/DBC;   %
                                cost = (DPCs-DBC);   %
                                %                                 cost = DPCs/DBC^2;
                                
                                tab_errors(count,:) = [pp trn_centers' ell1 ell2 DPCs DBC cost];
                                clear id1 id2 DPCs DBC cost K
                                count = count + 1;
                            end
                        end
                        clear count
                end
                %%%%%
                if strcmp(OO.method,'kernelDiff')
                    qmin = find(min(tab_errors(:,end)) == tab_errors(:,end)); % (tab_errors(:,2) == tab_errors(:,3))??
                    chosen_param = tab_errors(qmin(1),1:4);
                    cent_ind = tab_errors(qmin,5:6);
                    qmax = find(max(tab_errors(:,end)) == tab_errors(:,end));
                    mm = tab_errors(qmax(1),1:4);
                    DBC = tab_errors(qmin(1),end-1);
                    K = Kdiff(trainSet,trainSet,chosen_param,OO.ktype);
                    [lab_temp, cent_temp, DPC_T, DBC_T] = knkmeans(K,OO.Nclust);
                    clear qmin qmax
                    if strcmp(OO.plots,'on')
                        %%%
%                         figure; scatter3(trainSet(:,4),trainSet(:,3),trainSet(:,2),30,lab_temp,'f');drawnow, pause
                        %%%
                        % Im_error = reshape(tab_errors(:,end),length(unique(tab_errors(:,1))),length(unique(tab_errors(:,3))));
                        Im_error = reshape(standardimg(tab_errors(:,end),'zeroone'),length(pars2),length(pars1));
                        % Im_error2 = reshape(tab_errors(:,7),length(pars2),length(pars1));
                        % Im_error3 = reshape(tab_errors(:,8),length(pars2),length(pars1));
                        Im_error4 = reshape(tab_errors(:,9),length(pars2),length(pars1));
                        Im_error5 = reshape(tab_errors(:,10),length(pars2),length(pars1));
                        % contourf(pars1,pars2,Im_error);
                        % figure; imagesc(Im_error2); axis image, axis xy
                        %  figure; imagesc(Im_error3); axis image, axis xy
                        figure; imagesc(Im_error4); axis image, axis xy
                        figure; imagesc(Im_error5); axis image, axis xy
                        figure; imagesc(Im_error); axis image, axis xy
                        set(gca,'YDir','normal'); axis image
                        title(sprintf('par %f %f %f %f',chosen_param));
                        hold on
                        [myerr mxerr] = find(min(min(Im_error)) == Im_error);
                        scatter(mxerr,myerr,200,'white','filled')
                        [Myerr Mxerr] = find(max(max(Im_error)) == Im_error);
                        scatter(Mxerr,Myerr,200,'green','filled')
                        % legend('min', 'max')
                    end
                else
                    qmin = find(min(tab_errors(:,end)) == tab_errors(:,end)); % (tab_errors(:,2) == tab_errors(:,3))??
                    chosen_param = tab_errors(qmin(1),1);
                    cent_ind = tab_errors(qmin(1),2:3);
                    qmax = find(max(tab_errors(:,end)) == tab_errors(:,end));
                    mm = tab_errors(qmax(1),1);
                    DBC = tab_errors(qmin(1),end-1);
                    clear qmin qmax
                    K = kernelmatrix(OO.ktype,trainSet',trainSet',chosen_param);
                    [lab_temp, cent_temp, DPC_T, DBC_T] = knkmeans(K,OO.Nclust);
                    if strcmp(OO.plots,'on')
                        %%%
%                         K = Kdiff(trainSet,trainSet,chosen_param,OO.ktype);
%                         [lab_temp, cent_temp, DPC_T, DBC_T] = knkmeans(K,OO.Nclust);
                        figure; scatter3(trainSet(:,4),trainSet(:,3),trainSet(:,2),30,lab_trn,'f');drawnow, pause
                        %%%
                        figure;
                        plot(tab_errors(:,1),tab_errors(:,end))
                        %                         axis([0 size(tab_errors,1) 0 3])
                        title(sprintf('par %f',chosen_param));
                        hold on
                        plot(tab_errors(:,1),tab_errors(:,end-1),'red')
                        plot(tab_errors(:,1),tab_errors(:,end-2),'black')
                        %                         axis([0 10 0 1.5])
                        % plot(tab_errors(:,1),[0; diff(tab_errors(:,end))],'green')
                        % plot(tab_errors(:,1),[0; 0; diff(diff(tab_errors(:,end)))],'.-g')
                        % legend('COST','DBC','DPC norm sum','approx 1st der','approx 2nd der')
                    end
                end
                
                clear mm h mxerr myerr Mxerr Myerr p1 p2 pp
                %%%%%% KERNEL BLOCK
                % K = (chosen_param)*kernelmatrix('RBF',trainSet(:,1:4)',trainSet(:,1:4)',prior1) + ...
                %     (1-chosen_param)*kernelmatrix('RBF',trainSet(:,5:8)',trainSet(:,5:8)',prior2);
                %%%%%%
                %     if strcmp(OO.method,'kernelDiff')
                %         K = Kdiff(trainSet,trainSet,chosen_param,OO.ktype);
                %         %         K = normalise(K);
                %     else
                %         K = kernelmatrix(OO.ktype,trainSet',trainSet',chosen_param);
                %     end
                %     [lab, cent_ind, DPC, DBC] = knkmeans(K,OO.Nclust); % son da dividere per il numero di pixel per cluster
                %%%%%% TEST BLOCK: cluster with chosen parameter
                
                % % % INIT
                sizeBlocks = 200;
                index = 1;
                centers = trainSet(cent_ind',:);
                nblocks = ceil(OO.ColSize(1)/sizeBlocks);
                dist = zeros(OO.ColSize(1),OO.Nclust);
                stopping = min(index+sizeBlocks-1,OO.ColSize(1));
                for j = 1:nblocks
                    data_part = WholeSet(index:stopping,:);
                    if strcmp(OO.method,'kernelDiff')
                        self = diag(Kdiff(data_part,data_part,chosen_param,OO.ktype));
                    else
                        self = diag(kernelmatrix(OO.ktype,data_part',data_part',chosen_param));
                    end
                    if mod(j,100) == 0
                        fprintf('block %i of %i\n',j,nblocks);
                    end
                    for zz=1:OO.Nclust
                        if strcmp(OO.method,'kernelDiff')
                            fraPunti = 1/(sum(lab_temp == zz)^2)*sum(sum(Kdiff(trainSet(lab_temp == zz,:),trainSet(lab_temp == zz,:),chosen_param,OO.ktype)));
                            Kji = sum(Kdiff(trainSet(lab_temp == zz,:),data_part,chosen_param,OO.ktype),1)';
                            dist(index:stopping,zz) = repmat(fraPunti,size(data_part,1),1) + self - 2/(sum(lab_temp == zz))*Kji;
                        else
                            fraPunti = 1/(sum(lab_temp == zz)^2)*sum(sum(kernelmatrix(OO.ktype,trainSet(lab_temp == zz,:)',trainSet(lab_temp == zz,:)',chosen_param)));
                            dist(index:stopping,zz) = repmat(fraPunti,size(data_part,1),1) + self - ...
                                2/(sum(lab_temp == zz))*sum(kernelmatrix(OO.ktype,trainSet(lab_temp == zz,:)',data_part',chosen_param))';
                        end
                    end
                    index = index+sizeBlocks;
                    stopping = min(index+sizeBlocks-1,size(WholeSet,1));
                end
                % dist = sqrt(dist); % solo per correttezza, ma non cambia l'ordine delle cose.
                [min_dist,lab_out]=min(dist,[],2);
                
                %%% QUESTO PER L'IDEA ONE CLASS
                %                 if max(cent_ind) > size(trainSet,1)/2
                %                     labels = dist(:,2) > DBC;
                %                 else
                %                     labels = dist(:,1) > DBC/2;
                %                 end
                %                 [temp,labels]=min(dist,[],2);
                %%%%% clean up
                clear temp index stopping fraPunti data_part nblocks sizeBlocks self
                
                if strcmp(OO.plots,'on')
                    figure; imshow(reshape(lab_out,OO.OriginalSize(1),OO.OriginalSize(2)),[])
                    figure; imshow(reshape(dist(:,1),OO.OriginalSize(1),OO.OriginalSize(2)),[]); title('dist cluster 1')
                    figure; imshow(reshape(dist(:,2),OO.OriginalSize(1),OO.OriginalSize(2)),[]); title('dist cluster 2')
                    tilefigs
                end
                sigma_prior = OO.kpar;
                eval(sprintf('trn%i = trainSet;',zz))
% figure; plot(lab_out); drawnow;tilefigs
% pause
                mkdir([DIR ID])
                save([DIR ID 'EXP_' OO.method '_ktype_' OO.ktypedisp '_' num2str(OO.NtrnPix) 'trn' num2str(i) '.mat'],...
                    'chosen_param','lab_out','dist','tab_errors','sigma_prior');
                
                %                 eval(sprintf('save([DIR ID ''Res_SK_'' OO.method ''%i.mat''],''chosen_param'',''labels'',''dist'',''tab_errors'',''sigma_prior'')',i))
            end
            Read_exp
            %                         figure; imagesc(reshape(labels,OO.OriginalSize(1),OO.OriginalSize(2))); axis image, axis off
%             tilefigs
            clear i j sigma_prioir s pars
        end
    end
end

% for i = 1:5
%     figure; imagesc(Map_sum(:,:,i)); axis image
% end
% tilefigs
% load('./imm/Greece_island_ROInew.mat')
% figure; imagesc(ROI_n); axis image
% clear
Compute_skill_scores
Create_accuracy_tabs


% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%