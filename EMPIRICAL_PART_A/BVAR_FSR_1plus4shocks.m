% BVAR_FSR_1plus4shocks.m
% Code for fast estimation of structural restrictions in Bayesian SVARs. Compared to the baseline model the same number
% of shocks is estimated but we only place sign restrictions on the financial shock. The first four shocks are not
% identified as macro shocks anymore, and as there are no sign restrictions these now become nuisance shocks.
%
% This code replicates panel (c) of Figure 3 of Section 4.1 of the paper Korobilis, Dimitris (2022), "A new algorithm 
% for structural restrictions in vector autoregressions". https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3557911, 
% forthcoming in European Economic Review
%----------------------------------------------------------------------------------------------------------------------- 
% Written by Dimitris Korobilis
% University of Glasgow
% This version: 07 July 2022
%-----------------------------------------------------------------------------------------------------------------------

close all;
clear;
clc;

% Add path of functions
addpath('data');
addpath('functions');

%-------------------------------PRELIMINARIES--------------------------------------
ngibbs     = 500000;            % Gibbs sampler iterations
nburn      = 0.1*ngibbs;      % Iterations to discard
nthin      = 100;              % Thinning factor

% Please choose:
p          = 5;         % p is number of lags in the VAR part                         

% Impulse response analysis
nhor       = 36;
%----------------------------- END OF PRELIMINARIES --------------------------------
tic;
%----------------------------------LOAD DATA----------------------------------------   
% variables are in the following order:
% 1: Adjusted TFP
% 2: Stock Prices
% 3: Consumption
% 4: Real Interest Rate
% 5: Hours Worked
data  = xlsread('database.xlsx');
Y     = data(:,1:6);
names = {'GDP';'Prices';'Interest Rate';'Investment/Output';'Stock Prices';'External Finance Premium'};
% [Y,mY,stdY]     = zscore(Y);

% ===================================| VAR EQUATION |==============================
[y,x,M,T,KK,K] = prepare_BVAR_matrices(Y,p);

%==================| 1. OLS estimates:
[beta_OLS,beta_OLS_vec,SIGMA_OLS] = getOLS(y,x,M,p);

%==================| 2. PRIORS:
% Stuff for impulse response analysis: matrix of sign restrictions
%1=positive, -1=negative, 0=zero restriction,  NaN=no restriction
      %S1   S2   S3   S4   S5
sg = [NaN  NaN  NaN  NaN    1;    % VAR EQ1
      NaN  NaN  NaN  NaN    1;    % VAR EQ2
      NaN  NaN  NaN  NaN    1;    % VAR EQ3
      NaN  NaN  NaN  NaN    1;    % VAR EQ4
      NaN  NaN  NaN  NaN    1;    % VAR EQ5
      NaN  NaN  NaN  NaN  NaN];   % VAR EQ6
       
sg0 = sg; 
sg0(isnan(sg0)) = 0;
shock_names = {'S1';'S2';'S3';'S4';'Financial'};
nor = numel(find(sg~=0));  % number of restrictions 
N   = length(shock_names); % number of structural shocks

% Priors on model parameters
if size(x,1)>size(x,2)
    est_alg = 2;
else
    est_alg = 1;
end

% Horseshoe prior
psi = cell(M,1); tau = zeros(M,1); sigma_sq = zeros(M,1);
for ieq = 1:M
    psi{ieq,1} = ones(KK,1);
    tau(ieq,1) = 1;
    sigma_sq(ieq,1) = 1;
end

lam_prmean  = sg0;
lam_prvar   = 100*ones(M,N);
lam_iprvar  = 1./lam_prvar;
index_vec = 1:N;

% Extract factors
yhat = y - x*beta_OLS;
[F,L] = extract(zscore(yhat),N);
F = F/chol(cov(F)) - mean(F/chol(cov(F)));

%========|STORAGE MATRICES:
SIGMA   = SIGMA_OLS;
iSIGMA  = inv(SIGMA);
choliSIGMA = chol(iSIGMA);
V = kron(iSIGMA,speye(T));
index_kron = find(V~=0);
V = speye(M*T);
index_diag = find(V~=0);
R = eye(M);
iR = ones(M,1);
beta_mat = zeros(KK,M);

beta_save  = zeros(ngibbs/nthin,KK,M);
SIGMA_save = zeros(ngibbs/nthin,M,M);
L_save     = zeros(ngibbs/nthin,M,N);
F_save     = zeros(ngibbs/nthin,T,N);
R_save     = zeros(ngibbs/nthin,M);

irf_sign   = zeros(ngibbs/nthin,M,N,nhor);
D_theta    = zeros(ngibbs/nthin,1);
%======================= BEGIN MCMC ESTIMATION =======================
tic;
counter = 0;
for iter = 1:ngibbs + nburn
    if mod(iter,1000)==0   
        disp(['This is iteration ' num2str(iter)])         
        disp([num2str(100*(iter/(ngibbs+nburn))) '% completed'])
        toc;
    end
    
    y_til = y - F*L';    
    % STEP 1: update BETA
    for ieq = 1:M
       [beta_i,psi{ieq,1},tau(ieq)] = horseshoe(y_til(:,ieq),x,psi{ieq,1},tau(ieq),R(ieq,ieq),est_alg,ieq);
       beta_mat(:,ieq) = beta_i;
    end
    B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
    rej_count = 0;
    while max(abs(eig(B))) > 0.999
        rej_count = rej_count + 1;
        if rej_count > 500; y_til  = y; end
        for ieq = 1:M
            [beta_i,psi{ieq,1},tau(ieq)] = horseshoe(y_til(:,ieq),x,psi{ieq,1},tau(ieq),R(ieq,ieq),est_alg,ieq);
            beta_mat(:,ieq) = beta_i;         
        end
        B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
    end
    
    % STEP 2: update SIGMA 
    yhat = y - x*beta_mat;  % Now yhat is the VAR residual that follows the factor model          
    
    % Sample shocks F
    Lmsg = bsxfun(@times,L,iR);   
    Veta1 = eye(N) + Lmsg'*L;
    temp = cholcov(Veta1); [Q,C] = qr(temp);   
    S = inv(C); Veta = S*S';                     % Veta = inv(Veta1)
    F_mean = yhat*Lmsg*Veta;                     % T x N    
    F = F_mean + normrnd(0,1,[T,N])*S';          % sample F in a block
    
    % Normalize factors/shocks to be N(0,I) structural shocks
    F = F-mean(F);
    F = F/chol(cov(F));
    
    % I am doing now a very inefficient element-by-element update of the posterior of lambda,
    % so I can make a straightfoward use of the univariate truncated normal generator 
    for i = 1:M
        Lvar = diag(lam_iprvar(i,:)) + (F'*F).*iR(i); % inverse of posterior covariance matrix of L(i,:)
        L_bar(i,:) = Lvar\(diag(lam_iprvar(i,:))*lam_prmean(i,:)' + (F'*yhat(:,i)).*iR(i));
        for j = 1:N
            whole_vec = Lvar(j,:)'.*(L(i,:) - L_bar(i,:));
            Lpostvar = 1./Lvar(j,j); ss = sqrt(Lpostvar);
            Lpostmean = L_bar(i,j) - sum(whole_vec(find(index_vec~=j)))*Lpostvar;            
            if sg(i,j)  == 1
                lower  = (0 - Lpostmean)/ss;   
                upper  = (10000 - Lpostmean)/ss;
                L(i,j) = trandn(lower,upper);
                L(i,j) = Lpostmean + ss*L(i,j);
            elseif sg(i,j)  == -1
                lower  = (-10000 - Lpostmean)/ss;
                upper  = (0 - Lpostmean)/ss;
                L(i,j) = trandn(lower,upper);
                L(i,j) = Lpostmean + ss*L(i,j);
            elseif sg(i,j) == 0
                L(i,j)  = 0;
            elseif isnan(sg(i,j))
                L(i,j)  = Lpostmean + ss*randn;
            end
        end
    end
    
    % draw R
    sse2    = sum((yhat - F*L').^2);
    R_1     = (1 + T);    
    R_2     = (.01 + sse2);
    iR      = gamrnd(R_1./2,2./R_2)';
    R       = diag(1./iR);
    
    SIGMA  = L*cov(F)*L' + R;

    if iter>nburn && mod(iter,nthin) == 0
        counter = counter + 1;
        beta_save(counter,:,:)  = beta_mat; 
        SIGMA_save(counter,:,:) = SIGMA;
        L_save(counter,:,:)     = L;
        F_save(counter,:,:)     = F;
        R_save(counter,:)       = diag(R);
        
        % ------------------| IMPULSE RESPONSE ANALYSIS |------------------
        % 1) Obtain propagation coefficients in companion form       
        B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];           
        shock = L*chol(cov(F));
                       
        % 2) Obtain IRFs
        IRF = zeros(M,N,nhor);
        for ihor = 1:nhor
            Blarge = B^(ihor-1);
            IRF(:,:,ihor) = Blarge(1:M,1:M)*shock;
        end
        
        % save draws
        irf_sign(counter,:,:,:) = IRF;

        % Calculate DIC        
        c = -T*M/2*log(2*pi);
        u = y - x*beta_mat - F*L'; u = reshape(u',T*M,1);
        D_theta(counter,:) = c - T/2*log(det(R)) - .5*u'*kron(speye(T),diag(iR))*u;         
    end
end
%======================== END MCMC ESTIMATION ========================

BB = squeeze(median(beta_save,1));
LL = squeeze(median(L_save,1));
RR = squeeze(median(R_save,1));  iRR = diag(1./RR);
FF = squeeze(median(F_save,1));

if N==1
    u = y - x*BB - FF'*LL;
else
    u = y - x*BB - FF*LL'; 
end
u = reshape(u',T*M,1);
D_theta_bar = c - T/2*log(det(diag(RR))) - .5*u'*kron(speye(T),iRR)*u;

DIC = -4*mean(D_theta) + 2*D_theta_bar;

% ------------------| SAVE RESULTS |---------------------
% Save results in .mat file
modname = 'VAR_4plus1shocks';
save(sprintf('%s_%g_%g_%g.mat',modname,p,M,N),'-mat');

% irf_sign = irf_sign.*repmat(stdY,counter,1,N,nhor);

IRFlow = squeeze(quantile(irf_sign,.16,1));
IRFhigh = squeeze(quantile(irf_sign,.84,1));
IRFmed = squeeze(quantile(irf_sign,.5,1));


% First plot IRFs from Furlanetto, Ravazzolo, Sarferaz (forthcoming)
load('impresp.mat');

FRSlow = squeeze(quantile(response.f,.16,4));
FRShigh = squeeze(quantile(response.f,.84,4));
FRSmed = squeeze(quantile(response.f,.50,4));

% Plot IRFs
ylimits = zeros(5,6,2);
for j = 1:5
    fullscreen = get(0,'ScreenSize');
    figure('Position',[0 0 fullscreen(3) fullscreen(4)]);  
    for i = 1:6
        subplot(3,2,i)
        shadedplot((1:nhor)',squeeze(FRSlow(i,j,:))',squeeze(FRShigh(i,j,:))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(FRSmed(i,j,:))','g','LineWidth',3);  
        xlim([1 nhor]); grid on;   
        ax = gca;
        ax.FontSize = 16;
        title(names(i));
        %xlabel(['Quarters']);
        ylimits(j,i,:) = ylim;
    end
    if j ~= 4
        suptitle(['Responses to a ' str2mat(shock_names(j)) ' shock']);
    else
        suptitle(['Responses to an ' str2mat(shock_names(j)) ' shock']);
    end
end


% Next plot model-based IRFs
for j = 1:5
    fullscreen = get(0,'ScreenSize');
    figure('Position',[0 0 fullscreen(3) fullscreen(4)]);  
    for i = 1:6
        subplot(3,2,i)
        shadedplot((1:nhor)',squeeze(IRFlow(i,j,1:nhor))',squeeze(IRFhigh(i,j,1:nhor))',[0.8 0.8 0.8],[0.8 0.8 0.8]) 
        hold all
        plot((1:nhor)',zeros(1,nhor),'r',(1:nhor)',squeeze(IRFmed(i,j,1:nhor))','g','LineWidth',3);  
        xlim([1 nhor]); grid on;   
        ax = gca;
        ax.FontSize = 16;
        title(names(i));
        %xlabel(['Quarters']);
        ylim(ylimits(j,i,:));
    end
    if j ~= 4
        suptitle(['Responses to a ' str2mat(shock_names(j)) ' shock']);
    else
        suptitle(['Responses to an ' str2mat(shock_names(j)) ' shock']);
    end
    pngname = ['FACTOR_1plus4' num2str(j) '.png'];
    saveas(gcf,pngname);
end


clc; toc;

