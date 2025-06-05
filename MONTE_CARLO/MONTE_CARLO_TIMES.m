% MONTE_CARLO_TIMES.m  
% Monte Carlo exercise that simulates from a VAR with factor structure on the disturbances. Unlike the main Monte Carlo 
% exercise, parameters in this exercise are randomly generated and for simplicity I assume p=1 lags in the VARs and only
% 60,000 MCMC iterations. Obviously computing times are indicative and will vary widely when changing the number of lags,
% number of iterations, or number of endogenous variables in the VAR. Also times will vary widely based on the machine 
% (PC) you are using. For that reason this code comes "as is", with no guarantees in performance.
%
%----------------------------------------------------------------------------------------------------------------------- 
% This code replicates Table 3 of Section 3.2 of the paper Korobilis, Dimitris (2022), "A new algorithm for structural 
% restrictions in vector autoregressions". https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3557911, forthcoming
% in European Economic Review
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
nMC        = 500;              % Number of Monte Carlo iterations

ngibbs     = 50000;            % Gibbs sampler iterations
nburn      = 0.2*ngibbs;       % Iterations to discard
nthin      = 50;               % Thinning factor                   

% Impulse response analysis
nhor       = 60;
%----------------------------- END OF PRELIMINARIES --------------------------------
tic;
%----------------------------------GENERATE DATA----------------------------------------   
DGP = 2;
TIMES = zeros(nMC,1);
index = 0;
for T = [200,500]
    for M = 15
        for k = [3,10]
            index = index+1;
            for irep = 1:nMC
                disp('*****************************************');
                disp(['This is MC iteration ' num2str(irep)]);
                disp('*****************************************');
    
                params.T = T; params.M = M; params.k = k; params.p = 1;
                params.B = 0.8*eye(params.M); params.L = 2*rand(params.M,params.k) - 1; params.R = 0.1*diag(rand(params.M,1));
                [Y, beta_sim, sigma_sim]   =  simSVARDGP(params);

                sg  = sign(params.L);
     
                tic;
                [beta_save1,SIGMA_save1,L_save1,F_save1,irf_sign1,DIC1a,DIC1b] = Gibbs_FSR_TIMES(Y,1,sg,nhor,ngibbs,nburn,nthin);
                times = toc;
    
                TIMES(irep,index)= times;
            end
        end
    end
end
modname = 'TIMES';
save(sprintf('%s_%g.mat',modname,params.M),'-mat');
clc
