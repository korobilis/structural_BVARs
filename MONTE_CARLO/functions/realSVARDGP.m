function [y, BETA, SIGMA, L, F, R, IRF]   =  realSVARDGP(DGP)

% Simulate from VAR with factor structure on the covariance matrix
% =================================================================
% The model is of the form
%  
%   y[t] = B[1] x y[t-1] + ... + B[p] x y[t-p] + v[t],
%   v[t] = L x F[t] + e[t] 
%
% where y[t] is n x 1, F[t] is k x 1 and L = n x k.
% =================================================================

% Load real US data
addpath('data')
A=xlsread('data_for_MC.xlsx','data','B2:O518');
tcode = A(1,:);
data = A(2:end,:);
for i = 1:size(data,2); Y(:,i) = transx(data(:,i),tcode(i)); Y(:,i) = adjout(Y(:,i),4.5,4); end  % Convert to rates
Y = Y(2:end,:);   % Correct number of observations after converting to rates
Y = zscore(Y);    % Standardize for stability

% Estimate VAR with factor structure on the residuals using OLS
p = 12; % lags
k = 3;  % Number of shocks/factors
[y,x,n,T,~,] = prepare_BVAR_matrices(Y,p);  % Create VAR data matrices for estimation
[B,~,~] = getOLS(y,x,n,p);           % Obtain OLS estimate of AR matrix
yhat = y - x*B;                      % OLS residual
[F,L] = extract(zscore(yhat),k);     % Factor model on OLS residuals
%L = L/diag(diag(L));                 % Standardize shocks to be one for each driving variable
R = diag(mean((yhat - F*L').^2));    % Covariance of idiosyncratic shocks

if DGP == 1
    S = L*L' + R;
elseif DGP == 2           
    F = randn(T+1000,k);
    % Normalize factors to be N(0,I) structural shocks
    F = (grams(F')');   % First orthogonalize
    F = F./std(F);      % Then standardize
    e_t = F*L' + randn(T+1000,n)*sqrt(R);
    S = cov(e_t);
end

% Assign VAR coefficients for simulation
BETA  = B;
SIGMA = S;

% Extract IRFs based on these coefficients
nhor = 60;
B = [BETA(2:end,:)'; eye(n*(p-1)) , zeros(n*(p-1),n)];
shock = L;                     
IRF = zeros(n,k,nhor);
for ihor = 1:nhor
    Blarge = B^(ihor-1);
    IRF(:,:,ihor) = Blarge(1:n,1:n)*shock;
end


%----------------------GENERATE--------------------------
% Set storage in memory for y
y =[0*rand(p,n) ; zeros(T+1000-p,n)];

% Now generate Y from VAR (params.p,PHI,PSI)
if DGP == 1
    for t = p:T+1000
        u = chol(SIGMA)'*randn(n,1);
        ylag = mlag2(y,p);
        y(t,:) = [1 ylag(t,:)]*BETA + u';
    end    
elseif DGP == 2
    for t = p:T+1000
        ylag = mlag2(y,p);
        y(t,:) = [1 ylag(t,:)]*BETA + e_t(t,:);
    end    
end
y = y(end-T+1:end,:);

