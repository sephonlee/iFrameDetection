function [reg_rupt,reg_val,sig_rupt,sig_val]=Regsig_gener(t,opt_regsig)
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
%%%%%%%%%%%%%%%%%%
% Generates a regression function s and a noise level function sigma, randomly or not
%
% Input
% t	: vector of observation instants
% opt_regsig : integer between 1 and 999
%
% 	the random frameworks (A,B,C) of Section 2 the supplementary material correspond to opt_regsig=1,2,3 respectively
%
% 	the deterministic frameworks correspond to three-digits numbers
% 	the first digit is for the regression function:
%		1 for s_1
%		2 for s_2
%		3 for s_3
%		4 for s_4 (=s_1 with jumps shifted from 1e-4 to the right, see Section 4 of the Supplementary material)
% 	the second digit is for the shape of the noise-level function:
%		1 for constant (sigma_c)
%		2 for piecewise-constant (sigma_{pc})
%		3 for sinusoidal (sigma_s)
% 	the third digit is for the value of the noise-level (for sigma_{pc} only)
%		1 for constant (sigma_c)
%		2 for piecewise-constant (sigma_{pc})
%		3 for sinusoidal (sigma_s)
% 	Examples: 	opt_regsig=223 means (s_2,sigma_{pc,3})
%			opt_regsig=411 means (s_4,sigma_c)

%
% Output
% reg_rupt : vector of the breakpoints of s: r_1 = 0 < r_2 < ... < r_K=1 
% reg_val : vecteur of the values of s on each of these intervals: the i-th value is the value of s over [r_i,r_{i+1}[, except for i=1 (value on (-\infty,r_1)) and i=K-1 (value on [r_{K-1},+\infty)) 
% sig_rupt : vector of the breakpoints of sigma (same convention as for s)
% sig_val : vecteur of the values of sigma on each of these intervals (same convention as for s)
%
% Note: the breakpoints do not necessarily coincide with the design points given by t
%
n=numel(t);


if (opt_regsig<100)%%% Cadre aleatoire
    switch opt_regsig
    case {1}
        %%% Parametres "fixes" pour ce cadre de simulation
        %% Pour la fonction de regression
        nb_rupt_min_reg=3;% nombre minimal de ruptures pour reg
        nb_rupt_max_reg=floor(sqrt(n));% nombre maximal de ruptures pour reg
        taille_min_mcx_abs_reg=5/n;% borne inf sur la taille d'un morceau pour reg (modulo le nombre de ruptures choisi)
        saut_min_reg=0.1;% en valeur absolue, amplitude minimale d'un saut pour reg
        saut_max_reg=1;% en valeur absolue, amplitude maximale d'un saut pour reg
        %
        %% Pour le niveau de bruit
        nb_rupt_min_sig=5;% nombre minimal de ruptures pour sig
        nb_rupt_max_sig=floor(sqrt(n));% nombre maximal de ruptures pour sig
        taille_min_mcx_abs_sig=5/n;% borne inf sur la taille d'un morceau pour sig (modulo le nombre de ruptures choisi)
        sig_min=0.05;% valeur minimale pour sig 
        sig_max=0.5;% valeur maximale pour sig

        %%% generation de s
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_reg et nb_rupt_max_reg
        nb_rupt_reg=floor(nb_rupt_min_reg+max(1,nb_rupt_max_reg-nb_rupt_min_reg+1)*rand(1));
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_reg=min(taille_min_mcx_abs_reg,1/(1+nb_rupt_reg));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_reg+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_reg+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid uniformes sur [0;1]
        U_vect=rand(1,nb_rupt_reg+1);
        tmp=(1 - (nb_rupt_reg+1)*taille_min_mcx_reg) * U_vect/sum(U_vect);
        reg_rupt=[0 cumsum(taille_min_mcx_reg+tmp)];
        %
        %%% valeur de reg sur chacun des morceaux
        reg_val=cumsum((saut_min_reg+(saut_max_reg-saut_min_reg)*rand(1,nb_rupt_reg+1)).*sign(rand(1,nb_rupt_reg+1)-.5));

        %%% generation de sigma
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_sig et nb_rupt_max_sig
        nb_rupt_sig=floor(nb_rupt_min_sig+max(1,nb_rupt_max_sig-nb_rupt_min_sig+1)*rand(1));
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_sig=min(taille_min_mcx_abs_sig,1/(1+nb_rupt_sig));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_sig+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_sig+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid uniformes sur [0;1]
        U_vect=rand(1,nb_rupt_sig+1);
        tmp=(1 - (nb_rupt_sig+1)*taille_min_mcx_sig) * U_vect/sum(U_vect);
        sig_rupt=[0 cumsum(taille_min_mcx_sig+tmp)];
        %
        %%% valeur de sigma sur chacun des morceaux:
        %%  iid uniforme entre sig_min et sig_max
        sig_val=sig_min+(sig_max-sig_min).*rand(1,nb_rupt_sig+1);


    case {2}
        %%% Parametres "fixes" pour ce cadre de simulation
        %% Pour la fonction de regression
        nb_rupt_min_reg=3;% nombre minimal de ruptures pour reg
        nb_rupt_max_reg=floor(sqrt(n));% nombre maximal de ruptures pour reg
        taille_min_mcx_abs_reg=5/n;% borne inf sur la taille d'un morceau pour reg (modulo le nombre de ruptures choisi)
        saut_min_reg=0.1;% en valeur absolue, amplitude minimale d'un saut pour reg
        saut_max_reg=1;% en valeur absolue, amplitude maximale d'un saut pour reg
        %
        %% Pour le niveau de bruit
        nb_rupt_min_sig=5;% nombre minimal de ruptures pour sig
        nb_rupt_max_sig=floor(sqrt(n));% nombre maximal de ruptures pour sig
        taille_min_mcx_abs_sig=5/n;% borne inf sur la taille d'un morceau pour sig (modulo le nombre de ruptures choisi)
        sig_min=0.05;% valeur minimale pour sig 
        sig_max=0.5;% valeur maximale pour sig

        %%% generation de s
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_reg et nb_rupt_max_reg
        nb_rupt_reg=floor(nb_rupt_min_reg+max(1,nb_rupt_max_reg-nb_rupt_min_reg+1)*rand(1));
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_reg=min(taille_min_mcx_abs_reg,1/(1+nb_rupt_reg));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_reg+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_reg+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid, chacune egale a la valeur absolue d'une v.a. bimodale
        %%% (= N(0,1) avec proba 1/2 et N(10,1) avec proba 1/2)
        tmpA(1,:)=abs(randn(1,nb_rupt_reg+1));%premiere gaussienne
        tmpA(2,:)=abs(10+randn(1,nb_rupt_reg+1));%seconde gaussienne
        tmpB=floor(1+2*rand(1,nb_rupt_reg+1));%variable pour choisir entre les deux
        U_vect=diag(tmpA( tmpB,:))';

        tmp=(1 - (nb_rupt_reg+1)*taille_min_mcx_reg) * U_vect/sum(U_vect);
        reg_rupt=[0 cumsum(taille_min_mcx_reg+tmp)];
        %
        %%% valeur de reg sur chacun des morceaux
        reg_val=cumsum((saut_min_reg+(saut_max_reg-saut_min_reg)*rand(1,nb_rupt_reg+1)).*sign(rand(1,nb_rupt_reg+1)-.5));

        %%% generation de sigma
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_sig et nb_rupt_max_sig
        nb_rupt_sig=floor(nb_rupt_min_sig+max(1,nb_rupt_max_sig-nb_rupt_min_sig+1)*rand(1));
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_sig=min(taille_min_mcx_abs_sig,1/(1+nb_rupt_sig));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_sig+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_sig+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid uniformes sur [0;1]
        U_vect=rand(1,nb_rupt_sig+1);
        tmp=(1 - (nb_rupt_sig+1)*taille_min_mcx_sig) * U_vect/sum(U_vect);
        sig_rupt=[0 cumsum(taille_min_mcx_sig+tmp)];
        %
        sig_val=sig_min+(sig_max-sig_min).*rand(1,nb_rupt_sig+1);

    case {3}
        %%%%%%%% L'intervalle est coupe en deux:
        %%%%%%%% 1ere partie: bruit faible, beaucoup de ruptures pour reg
        %%%%%%%% 2eme partie: bruit fort, peu de ruptures pour reg
        %%%%%%%%
        %%%%%%%% Les ruptures pour sigma sont prises nombreuses, independamment de reg
        %%%%%%%% Les ruptures pour reg sont choisies comme pour {7-10}, avec une distribution bimodale favorisant les plages courtes et les plages longues
        %%%%%%%%
        %%%%%%%%
        %%% Parametres "fixes" pour ce cadre de simulation
        %% Pour la coupure en deux regimes
        coupure=1/2;% point de coupure (deterministe) entre les 2 regimes
        %
        %
        %% Pour la fonction de regression
        %% Seconde partie (peu de ruptures pour reg)
        nb_rupt_min_reg_deux=0;% nombre minimal de ruptures pour reg
        nb_rupt_max_reg_deux=floor((floor(sqrt(n))-1)/3);% nombre maximal de ruptures pour reg
        %% Premiere partie (beaucoup de ruptures pour reg)
        nb_rupt_min_reg_un=2;% nombre minimal de ruptures pour reg
        nb_rupt_max_reg_un=floor(sqrt(n))-1-nb_rupt_max_reg_deux;% nombre maximal de ruptures pour reg
        %
        %% Les deux parties
        taille_min_mcx_abs_reg=5/n;% borne inf sur la taille d'un morceau pour reg (modulo le nombre de ruptures choisi)
        saut_min_reg=0.1;% en valeur absolue, amplitude minimale d'un saut pour reg
        saut_max_reg=1;% en valeur absolue, amplitude maximale d'un saut pour reg
        %
        %
        %% Pour le niveau de bruit
        nb_rupt_min_sig=floor(2*sqrt(n)/3);% nombre minimal de ruptures pour sig
        nb_rupt_max_sig=floor(4*sqrt(n)/3);% nombre maximal de ruptures pour sig
        taille_min_mcx_abs_sig=2/n;% borne inf sur la taille d'un morceau pour sig (modulo le nombre de ruptures choisi)
        sig_min_petit=0.025;% valeur minimale pour sig (1ere partie; bruit faible)
        sig_max_petit=0.2;% valeur maximale pour sig (1ere partie; bruit faible)
        sig_min_grand=0.1;% valeur minimale pour sig (2nde partie; bruit fort)
        sig_max_grand=0.8;% valeur maximale pour sig (2nde partie; bruit fort)

        %%% generation de s
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_reg et nb_rupt_max_reg (pour chaque partie)
        nb_rupt_reg_un=floor(nb_rupt_min_reg_un+max(1,nb_rupt_max_reg_un-nb_rupt_min_reg_un+1)*rand(1));
        nb_rupt_reg_deux=floor(nb_rupt_min_reg_deux+max(1,nb_rupt_max_reg_deux-nb_rupt_min_reg_deux+1)*rand(1));
        nb_rupt_reg=nb_rupt_reg_un+nb_rupt_reg_deux+1;
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_reg_un=min(taille_min_mcx_abs_reg,1/(1+nb_rupt_reg_un));
        taille_min_mcx_reg_deux=min(taille_min_mcx_abs_reg,1/(1+nb_rupt_reg_deux));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_reg+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_reg+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid, chacune egale a la valeur absolue d'une v.a. bimodale
        %%% (= N(0,1) avec proba 1/2 et N(0,10) avec proba 1/2)
        %% Premiere partie
        tmpA_un=zeros(2,nb_rupt_reg_un+1);
        tmpA_un(1,:)=abs(randn(1,nb_rupt_reg_un+1));%premiere gaussienne
        tmpA_un(2,:)=abs(10+randn(1,nb_rupt_reg_un+1));%seconde gaussienne
        tmpB_un=floor(1+2*rand(1,nb_rupt_reg_un+1));%variable pour choisir entre les deux
        U_vect_un=diag(tmpA_un( tmpB_un,:))';
        tmp_un = (1 - (nb_rupt_reg_un+1)*taille_min_mcx_reg_un) * U_vect_un/sum(U_vect_un);
        reg_rupt_un = coupure*[0 cumsum(taille_min_mcx_reg_un+tmp_un)];
        %
        %% Seconde partie
        tmpA_deux=zeros(2,nb_rupt_reg_deux+1);
        tmpA_deux(1,:)=abs(randn(1,nb_rupt_reg_deux+1));%premiere gaussienne
        tmpA_deux(2,:)=abs(10+randn(1,nb_rupt_reg_deux+1));%seconde gaussienne
        tmpB_deux=floor(1+2*rand(1,nb_rupt_reg_deux+1));%variable pour choisir entre les deux
        U_vect_deux=diag(tmpA_deux( tmpB_deux,:))';
        tmp_deux = (1 - (nb_rupt_reg_deux+1)*taille_min_mcx_reg_deux) * U_vect_deux/sum(U_vect_deux);
        reg_rupt_deux = coupure + (1-coupure)*cumsum(taille_min_mcx_reg_deux+tmp_deux);
        %
        %% Bilan
        reg_rupt=[reg_rupt_un reg_rupt_deux];
        %
        %%% valeur de reg sur chacun des morceaux
        reg_val=cumsum((saut_min_reg+(saut_max_reg-saut_min_reg)*rand(1,nb_rupt_reg+1)).*sign(rand(1,nb_rupt_reg+1)-.5));

        %%% generation de sigma
        %
        %%% nombre de ruptures: uniforme entre nb_rupt_min_sig et nb_rupt_max_sig
        nb_rupt_sig=floor(nb_rupt_min_sig+max(1,nb_rupt_max_sig-nb_rupt_min_sig+1)*rand(1));
        %%% taille minimale d'un morceau (deterministe)
        taille_min_mcx_sig=min(taille_min_mcx_abs_sig,1/(1+nb_rupt_sig));
        %
        %%% position des ruptures
        %%% On positionne les ruptures en calculant la taille des nb_rupt_sig+1 intervalles successifs
        %%% La taille du i-eme intervalle est taille_min_mcx + epsilon_i 
        %%% ou epsilon_i = (1 - (nb_rupt_sig+1)*taille_min_mcx) * U_i / (sum_j U_j)
        %%% et les U_j sont iid uniformes sur [0;1]
        U_vect=rand(1,nb_rupt_sig+1);
        tmp=(1 - (nb_rupt_sig+1)*taille_min_mcx_sig) * U_vect/sum(U_vect);
        sig_rupt=[0 cumsum(taille_min_mcx_sig+tmp)];
        %
        %%% valeur de sigma sur chacun des morceaux:
        %% on commence par calculer le vecteur de sig_min et sig_max (puisque sigma change)
        sig_min=sig_min_grand*ones(1,nb_rupt_sig+1);
        sig_max=sig_max_grand*ones(1,nb_rupt_sig+1);
        sig_min(sig_rupt(2:end)<coupure)=sig_min_petit;
        sig_max(sig_rupt(2:end)<coupure)=sig_max_petit;

        sig_val=sig_min+(sig_max-sig_min).*rand(1,nb_rupt_sig+1);

    otherwise
        error('Erreur : cadre aleatoire inexistant');
    end
else%%% Cadre deterministe
    reg_num=floor(opt_regsig/100);
    sig_typ=floor(opt_regsig/10)-10*reg_num;
    sig_num=floor(opt_regsig)-100*reg_num-10*sig_typ;
    switch reg_num
        case 1
            reg_rupt=(0:0.2:1);
            reg_val=[1 0 1 0 1];
        case 2
            reg_rupt=[0 , 0.35 , 0.55 , 0.7 , 0.8 , 1];
            reg_val=[0.5 , -0.5 , 0 , 0.25  ,0.5];
        case 3
            reg_rupt=[0 , 0.16 , 0.26 , 0.33 , 0.38 , 0.48 , 0.53 , 0.62 , 0.82 , 0.87 , 1];
            reg_val=[0.5 , -0.5 , 0 , 0.25 , 0.5 , 1 , 0.66 , 0.33 , 0.16 , 0];
        case 4
            reg_rupt=[0 0.2001 0.4001 0.6001 0.8001 1];
            reg_val=[1 0 1 0 1];
        otherwise
            error('Fonction de regression deterministe inconnue');
    end%switch reg_num
    %
    switch sig_typ
        case 1
            sig_rupt=[0 , 1];
            sig_val=0.25;
        case 2
            sig_rupt=[0 , 1/3 , 1];
            sig_val=[0.2 , 0.05];
            if sig_num==2
                sig_val=2*sig_val;
            end
            if sig_num==3
                sig_val=2.5*sig_val;
            end
        case 3
            sig_rupt=[t 2];
            sig_val=0.5*sin(t*pi/4);
        otherwise
            error('Fonction variance deterministe inconnue');
    end%switch sig_typ
            
end
    
