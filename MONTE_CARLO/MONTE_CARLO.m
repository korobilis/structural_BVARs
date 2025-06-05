% MONTE_CARLO.m  
% Monte Carlo exercise that simulates from a VAR with factor structure on the disturbances. Parameters for this VAR come
% from a real US dataset that is described in Appendix B.2 of the paper (see directory ~\data\data_for_MC.xlsx)
%
%----------------------------------------------------------------------------------------------------------------------- 
% This code replicates Figure 1 and Table 2 of Section 3.1, as well as Figures C1-C4 of the online supplement, of the 
% paper Korobilis, Dimitris (2022), "A new algorithm for structural restrictions in vector autoregressions". 
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3557911, forthcoming in European Economic Review
%
%-----------------------------------------------------------------------------------------------------------------------
% Written by Dimitris Korobilis
% University of Glasgow
% This version: 07 July 2022
%-----------------------------------------------------------------------------------------------------------------------

clear;
clc;

% Add path of data and functions
addpath('functions');
addpath('data');

%-------------------------------PRELIMINARIES--------------------------------------
nMC        = 500;               % Number of Monte Carlo iterations

ngibbs     = 500000;             % Gibbs sampler iterations
nburn      = 0.1*ngibbs;       % Iterations to discard
nthin      = 100;               % Thinning factor                   

% Impulse response analysis
nhor       = 60;
%----------------------------- END OF PRELIMINARIES --------------------------------
tic;
%----------------------------------GENERATE DATA----------------------------------------   
DGP = 2;
BETA   = zeros(169,14,5,nMC);
LAMBDA = zeros(14,3,5,nMC);
IRF    = zeros(14,3,nhor,5,nMC);
DIC    = zeros(5,nMC);
for irep = 1:nMC
    disp('*****************************************');
    disp(['This is MC iteration ' num2str(irep)]);
    disp('*****************************************');
    %[Y, beta_sim, sigma_sim, L_sim, F_sim, R_sim, IRF_sim]   =  realSVARDGP(DGP);
    [Y, beta_sim, sigma_sim, L_sim, F_sim, R_sim, IRF_sim]   =  realVARDGP;

    % matrix of sign restrictions
          %S1   S2   S3 
    sg = [  1,  -1,  -1;    % VAR EQ1
            1,   1,  -1;    % VAR EQ2
          NaN, NaN,   1;    % VAR EQ3
          NaN, NaN,  -1;    % VAR EQ4
          NaN, NaN, NaN;    % VAR EQ5
          NaN, NaN,  -1;    % VAR EQ6
          NaN, NaN,  -1;    % VAR EQ7
           -1,  -1,  -1;    % VAR EQ8
           -1,   1,   1;    % VAR EQ9
            1,  -1,  -1;    % VAR EQ10
            1,  -1,  -1;    % VAR EQ11
            1,   1,  -1;    % VAR EQ12
            1,   1,  -1;    % VAR EQ13
            1,   1,  -1];   % VAR EQ14
     sg = nan(14,4);
     sg0 = sg;
     sg0(isnan(sg0))=0;
     
     %% CASE 1: Correct specification of data, lags, restrictions, shocks
     [beta_save1,SIGMA_save1,L_save1,F_save1,irf_sign1,DIC1a,DIC1b] = Gibbs_FSR(Y,12,sg,nhor,ngibbs,nburn,nthin);

     %% CASE 2: Smaller VAR without extra data on output and prices
     [beta_save2,SIGMA_save2,L_save2,F_save2,irf_sign2,DIC2a,DIC2b] = Gibbs_FSR(Y(:,1:8),12,sg(1:8,:),nhor,ngibbs,nburn,nthin);

     %% CASE 3: Misspecification of lag structure
     [beta_save3,SIGMA_save3,L_save3,F_save3,irf_sign3,DIC3a,DIC3b] = Gibbs_FSR(Y,2,sg,nhor,ngibbs,nburn,nthin);
        
     %% CASE 4: Misspecification of shocks (estimate 1 less)
     [beta_save4,SIGMA_save4,L_save4,F_save4,irf_sign4,DIC4a,DIC4b] = Gibbs_FSR(Y,12,sg(:,1:2),nhor,ngibbs,nburn,nthin);
        
     %% CASE 5: Misspecification of shocks (estimate 1 more)
     sgsupp = [1;1;1;-1;1;NaN;NaN;-1;1;1;1;1;1;1]; sgsupp0 = sgsupp; sgsupp0(isnan(sgsupp0))=0;  % additional vector of restrictions for shock No4
     [beta_save5,SIGMA_save5,L_save5,F_save5,irf_sign5,DIC5a,DIC5b] = Gibbs_FSR(Y,12,[sg,sgsupp],nhor,ngibbs,nburn,nthin);

     %% SAVE MATRICES 
     BETA(:,:,1,irep)      = squeeze(median(beta_save1));
     BETA(1:97,1:8,2,irep) = squeeze(median(beta_save2));
     BETA(1:29,:,3,irep)   = squeeze(median(beta_save3));
     BETA(:,:,4,irep)      = squeeze(median(beta_save4));
     BETA(:,:,5,irep)      = squeeze(median(beta_save5));
        
     LAMBDA(:,:,1,irep)    = squeeze(median(L_save1));
     LAMBDA(1:8,:,2,irep)  = squeeze(median(L_save2));
     LAMBDA(:,:,3,irep)    = squeeze(median(L_save3));
     LAMBDA(:,1:2,4,irep)  = squeeze(median(L_save4));
     LAMBDA(:,:,5,irep)    = squeeze(median(L_save5(:,:,1:3)));
        
     IRF(:,:,:,1,irep)     = squeeze(median(irf_sign1));
     IRF(1:8,:,:,2,irep)   = squeeze(median(irf_sign2));
     IRF(:,:,:,3,irep)     = squeeze(median(irf_sign3));
     IRF(:,1:2,:,4,irep)   = squeeze(median(irf_sign4));
     IRF(:,:,:,5,irep)     = squeeze(median(irf_sign5(:,:,1:3,:)));
        
     DIC(1,irep) = DIC1a; DIC(2,irep) = DIC2a; DIC(3,irep) = DIC3a; DIC(4,irep) = DIC4a; DIC(5,irep) = DIC5a;
     DIC2(1,irep) = DIC1b; DIC2(2,irep) = DIC2b; DIC2(3,irep) = DIC3b; DIC2(4,irep) = DIC4b; DIC2(5,irep) = DIC5b;
end

save MONTE_CARLO.mat;
IRF(1:2,:,:,:,:) = cumsum(IRF(1:2,:,:,:,:),3);
IRF_sim(1:2,:,:) = cumsum(IRF_sim(1:2,:,:),3);


IRFlow = squeeze(quantile(IRF,.05,5));
IRFhigh = squeeze(quantile(IRF,.95,5));
IRFmed = squeeze(quantile(IRF,.5,5));


%% PLOT STUFF
fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for i = 1:3
    for j = 1:3
        subplot(3,3,(j-1)*3+i)
        shadedplot((1:nhor)',squeeze(IRFlow(j,i,:,1))',squeeze(IRFhigh(j,i,:,1))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(j,i,:,1))','g','LineWidth',3);
        hold all
        plot(squeeze(IRF_sim(j,i,:)),'LineWidth',3,'Color','k','LineStyle',':')    
        xlim([1 nhor]); grid on;
        if j==1
            title(['s_{' num2str(i) '}'],'FontWeight','bold');
        end
        if i==1       
            ylabel(['y_{' num2str(j) '}'],'FontWeight','bold');
        end
        ax = gca;
        ax.FontSize = 16;
    end
end
suptitle('CASE 1: Correct specification of data, lags, restrictions, shocks')
saveas(gcf,'CASE_1.png');


figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for i = 1:3
    for j = 1:3
        subplot(3,3,(j-1)*3+i)
        shadedplot((1:nhor)',squeeze(IRFlow(j,i,:,2))',squeeze(IRFhigh(j,i,:,2))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(j,i,:,2))','g','LineWidth',3);
        hold all
        plot(squeeze(IRF_sim(j,i,:)),'LineWidth',3,'Color','k','LineStyle',':')
        xlim([1 nhor]); grid on;
        if j==1
            title(['s_{' num2str(i) '}'],'FontWeight','bold');
        end
        if i==1       
            ylabel(['y_{' num2str(j) '}'],'FontWeight','bold');
        end
        ax = gca;
        ax.FontSize = 16;        
    end
end
suptitle('CASE 2: Smaller VAR without additional data on output and prices')
saveas(gcf,'CASE_2.png');

figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for i = 1:3
    for j = 1:3
        subplot(3,3,(j-1)*3+i)
        shadedplot((1:nhor)',squeeze(IRFlow(j,i,:,3))',squeeze(IRFhigh(j,i,:,3))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(j,i,:,3))','g','LineWidth',3);
        hold all
        plot(squeeze(IRF_sim(j,i,:)),'LineWidth',3,'Color','k','LineStyle',':')
        xlim([1 nhor]); grid on;
        if j==1
            title(['s_{' num2str(i) '}'],'FontWeight','bold');
        end
        if i==1       
            ylabel(['y_{' num2str(j) '}'],'FontWeight','bold');
        end
        ax = gca;
        ax.FontSize = 16;        
    end
end
suptitle('CASE 3: Misspecification of lag structure')
saveas(gcf,'CASE_3.png');

figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for i = 1:2
    for j = 1:3
        subplot(3,2,(j-1)*2+i)
        shadedplot((1:nhor)',squeeze(IRFlow(j,i,:,4))',squeeze(IRFhigh(j,i,:,4))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(j,i,:,4))','g','LineWidth',3);
        hold all
        plot(squeeze(IRF_sim(j,i,:)),'LineWidth',3,'Color','k','LineStyle',':')
        xlim([1 nhor]); grid on;
        if j==1
            title(['s_{' num2str(i) '}'],'FontWeight','bold');
        end
        if i==1       
            ylabel(['y_{' num2str(j) '}'],'FontWeight','bold');
        end
        ax = gca;
        ax.FontSize = 16;        
    end
end
suptitle('CASE 4: Misspecification of shocks (estimate 1 less)')
saveas(gcf,'CASE_4.png');

figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for i = 1:3
    for j = 1:3
        subplot(3,3,(j-1)*3+i)
        shadedplot((1:nhor)',squeeze(IRFlow(j,i,:,5))',squeeze(IRFhigh(j,i,:,5))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(j,i,:,5))','g','LineWidth',3);
        hold all
        plot(squeeze(IRF_sim(j,i,:)),'LineWidth',3,'Color','k','LineStyle',':')
        xlim([1 nhor]); grid on;
        if j==1
            title(['s_{' num2str(i) '}'],'FontWeight','bold');
        end
        if i==1       
            ylabel(['y_{' num2str(j) '}'],'FontWeight','bold');
        end
        ax = gca;
        ax.FontSize = 16;        
    end
end
suptitle('CASE 5: Misspecification of shocks (estimate 1 more)')
saveas(gcf,'CASE_5.png');


DD1 = squeeze(mean(DIC,2));
DD2 = squeeze(mean(DIC2,2));

disp('Table 1: DIC values based on 14 and 8 variables')
disp(DD1)
disp(DD2)

clc; toc;

