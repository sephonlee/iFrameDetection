%%%%%%%%%% Demo for the package 'CHPTCV'
%%%%%%%%%%
%%%%%%%%%% Illustrates the simulation study of the paper
%%%%%%%%%% 'Segmentation of the mean of heteroscedastic data via cross-validation'
%%%%%%%%%% by Sylvain Arlot and Alain Celisse (2010). Statistics and Computing, 1--20. DOI: 10.1007/s11222-010-9196-x. 
%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tunable Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Data type
%%% 1 : your own real data (uses 'nom_load')
%%% 2 : simulated data (uses 'options.opt_regsig', 'options.loi_err' and 'n')
z_datatype=2;
    nom_load='myfavoritedata.mat';
    options.opt_regsig=223; % see Regsig_gener.m for the correspondence with frameworks considered in the paper
    options.loi_err=1; % Distribution of the errors epsilon_i (1: Gaussian; 2: Exponential)
    n=100; % number of observations

%%% Initial plot of data
z_dataplot=1;% 1 : start by plotting data (+ regression function if available) 0: do nothing

%%% Procedures to be considered
%%% each of the variables z_XXX below correspond to one procedure
%%% z_XXX = 1 means the procedure will be launched on Data
%%% z_XXX = 0 means the procedure will not be launched
z_LOOVF=1;
z_ERMVF=1;
z_ERMBM=1;
z_ZS=1;
z_PML=1;

%%% Parameters of the procedures
infos.delta=1;% (Minimal number of design points between two breakpoints)-1
V=5;% How many folds for ERMVF and LooVF ?
Dimmax=floor(0.9*(n-n/V)/(1+infos.delta));% Maximal dimension considered for all procedures (default choice: our advice for LooVF or ERMVF)
infos.threshold=floor(0.75*Dimmax);% value of the dimension threshold used in ERMBMthr and PML for the calibration of the penalty 

%%% Plot crit_2(D) as a function of D
z_plotcrit2=1;% 1 : Plot crit_2(D) as a function of D for each procedure considered ; 0 : do nothing

%%% Visualize the partition selected for a given number of breakpoints (uses 'D_chosen')
z_visuDchosen=0;% 1 : Plot data + selected partition for each procedure considered ; 0 : do nothing
    D_chosen=5;% Dimension of the partition to be plotted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of tunable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%
%%%%%%%%%% 1. Data generation / loading
%%%%%%%%%%

switch z_datatype
    case 1

        %%%%%%%%%% 1.a. Analyze your own data set

        %%% load your data set

        load(nom_load);

        %%% if necessary, put your sample Y_1, ... , Y_n in the vector 'Data'
        %%% and the time values t_1, ... , t_n in the vector 't'
        %
        %Data=data(50:150);

        t=(1:numel(Data));
        %t=50:150;

        %%% for simulated data, put the values of the signal s(t_1), ... , s(t_n) in the vector 'Reg'
        %%% otherwise, just put zeros in Reg

        Reg=zeros(size(Data));

    case 2
        %%%%%%%%%% 1.b. Generate data according to one of the settings of the paper

        %%% Computes the corresponding pair (s,sigma) and put them in 'options'

        t=(1:n)/n;
        [options.reg_rupt,options.reg_val,options.sig_rupt,options.sig_val]=Regsig_gener(t,options.opt_regsig);

        %%% Generates data

        [Data Reg]=DataGener(t,options.reg_rupt,options.reg_val,options.sig_rupt,options.sig_val,options.loi_err);

    otherwise
        fprintf('Error: z_datatype\n');
end

n=numel(Data);

%%%%%%%%%%
%%%%%%%%%% 2. (Optional) Data plot
%%%%%%%%%%

if z_dataplot==1

    %%%%%%%%%% 2.a Plot the data set 
    
    figure(1);
    newplot;
    plot(t,Data,'b.');
    legend('data');
    xlabel('t');
    ylabel('Y');

    %%%%%%%%%% 2.b For simulated data, plot the regression function

    if z_datatype==2
        figure(1);
        set(gca,'NextPlot','add');
        plot((0:1e-5:1),eval_histo(options.reg_val,options.reg_rupt,(0:1e-5:1)),'r-');
        legend('data','regression function');
    end

end

%%%%%%%%%%
%%%%%%%%%% 3. Run the procedures
%%%%%%%%%%


%%%%%%%%%% 3.a LooVF

if z_LOOVF==1
    [D_LOOVF, mu_LOOVF, rupt_LOOVF, crit2VF_1LOO, rupt_1LOO, mu_1LOO]=proc_LOOVF(Data, V, Dimmax, infos);
end

%%%%%%%%%% 3.b ERMVF

if z_ERMVF==1
    [D_ERMVF, mu_ERMVF, rupt_ERMVF, crit2VF_1ERM, rupt_1ERM, mu_1ERM]=proc_ERMVF(Data, V, Dimmax, infos);
end
    

%%%%%%%%%% 3.c ERMBM

if z_ERMBM==1
    [D_ERMBMest, mu_ERMBMest, rupt_ERMBMest, crit2BMest_1ERM, D_ERMBMthr, mu_ERMBMthr, rupt_ERMBMthr, crit2BMthr_1ERM, D_ERMBMmax, mu_ERMBMmax, rupt_ERMBMmax, crit2BMmax_1ERM, rupt_1ERM, mu_1ERM]=proc_ERMBM(Data, Dimmax, infos);
end

%%%%%%%%%% 3.d ZS

if z_ZS==1
    [D_ZS, mu_ZS, rupt_ZS, crit2ZS, rupt_1ZS, mu_1ZS]=proc_ZS(Data, Dimmax, infos);
end

%%%%%%%%%% 3.e PML

if z_PML==1
    [D_PML, mu_PML, rupt_PML, crit2PML, rupt_1PML, mu_1PML]=proc_PML(Data, Dimmax, infos);
end



%%%%%%%%%%
%%%%%%%%%% 4. Visualize the estimators proposed by each procedure
%%%%%%%%%%


if z_LOOVF==1
    figure(11);newplot;
    plot(t,mu_LOOVF,'r-',t,Data,'b.');
    xlabel('t');ylabel('Y');
    legend('[LOO,VF]','data');
end

if z_ERMVF==1
    figure(12);newplot;
    plot(t,mu_ERMVF,'r-',t,Data,'b.');
    xlabel('t');ylabel('Y');
    legend('[ERM,VF]','data');
end

if z_ERMBM==1
    figure(13);newplot;
    plot(t,mu_ERMBMest,'r-',t,Data,'b.');
    hold on;
    plot(rupt_ERMBMest, mu_ERMBMest(rupt_ERMBMest), '^c')
    xlabel('t');ylabel('Y');
    legend('[ERM,BM]','data');
end

if z_ZS==1
    figure(14);newplot;
    plot(t,mu_ZS,'r-',t,Data,'b.');
    xlabel('t');ylabel('Y');
    legend('ZS','data');
end

if z_PML==1
    figure(15);newplot;
    plot(t,mu_PML,'r-',t,Data,'b.');
    xlabel('t');ylabel('Y');
    legend('PML','data');
end



%%%%%%%%%%
%%%%%%%%%% 5. Compute the performance of each procedure (only if the
%%%%%%%%%% true signal is known, i.e., for simulated data sets)
%%%%%%%%%%

if z_LOOVF==1
    risk_LOOVF=mean((mu_LOOVF-Reg).^2);
end

if z_ERMVF==1
    risk_ERMVF=mean((mu_ERMVF-Reg).^2);
end

if z_ERMBM==1
    risk_ERMBMest=mean((mu_ERMBMest-Reg).^2);
    risk_ERMBMthr=mean((mu_ERMBMthr-Reg).^2);
    risk_ERMBMmax=mean((mu_ERMBMmax-Reg).^2);
end

if z_ZS==1
    risk_ZS=mean((mu_ZS-Reg).^2);
end

if z_PML==1
    risk_PML=mean((mu_PML-Reg).^2);
end


%%%%%%%%%%
%%%%%%%%%% 6. Plot the value of the criterion used for choosing D
%%%%%%%%%%

if z_plotcrit2==1

    if z_LOOVF==1
        figure(21);newplot;
        plot(1:Dimmax,crit2VF_1LOO,'k-');
        xlabel('Dimension D');ylabel('crit_{2,VF}(D)');
        legend('[LOO,VF]');
    end

    if z_ERMVF==1
        figure(22);newplot;
        plot(1:Dimmax,crit2VF_1ERM,'k-');
        xlabel('Dimension D');ylabel('crit_{2,VF}(D)');
        legend('[ERM,VF]');
    end

    if z_ERMBM==1
        figure(23);newplot;
        plot(1:Dimmax,crit2BMest_1ERM,'r-',1:Dimmax,crit2BMthr_1ERM,'b--',1:Dimmax,crit2BMmax_1ERM,'g-.');
        xlabel('Dimension D');ylabel('crit_{2,VF}(D)');
        legend('[ERM,BM]','[ERM,BM thr.]','[ERM,BM max]');
    end

    if z_ZS==1
        figure(24);newplot;
        plot(1:Dimmax,-crit2ZS,'k-');
        xlabel('Dimension D');ylabel('-crit_{2,ZS}(D)');
        legend('ZS');
    end

    if z_PML==1
        figure(25);newplot;
        plot(1:Dimmax,crit2PML,'k-');
        xlabel('Dimension D');ylabel('crit_{2,PML}(D)');
        legend('PML');
    end
end


%%%%%%%%%%
%%%%%%%%%% 7. Visualize the partitions selected for various values of D inside each procedure
%%%%%%%%%%

if z_visuDchosen==1

    %%% Plots

    if z_LOOVF==1
        figure(31);newplot;
        plot(t,mu_1LOO(D_chosen,:),'r-',t,Data,'b.');
        xlabel('t');ylabel('Y');
        legend(['LOO with D=' num2str(D_chosen)],'data');
    end

    if (z_ERMVF==1 || z_ERMBM==1)
        figure(32);newplot;
        plot(t,mu_1ERM(D_chosen,:),'r-',t,Data,'b.');
        xlabel('t');ylabel('Y');
        legend(['ERM with D=' num2str(D_chosen)],'data');
    end


    if (z_PML==1)
        figure(33);newplot;
        plot(t,mu_1PML(D_chosen,:),'r-',t,Data,'b.');
        xlabel('t');ylabel('Y');
        legend(['PML with D=' num2str(D_chosen)],'data');
    end

end


%%%%%%%%%%
%%%%%%%%%% 
%%%%%%%%%%
