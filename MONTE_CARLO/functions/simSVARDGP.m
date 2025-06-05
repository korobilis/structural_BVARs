function [y, BETA, SIGMA]   =  simSVARDGP(params)

% Simulate from VAR with factor structure on the covariance matrix
% =================================================================
% The model is of the form
%  
%   y[t] = B x y[t-1] + v[t],
%   v[t] ~ N(0, S) 
%
% where y[t] is M x 1.
%
% INPUT: vector params with all the parameter values for the DGP
%    params.T    Number of observations
%    params.M    Number of endogenous variables
%    params.k    Number of structural shocks ("factors")
%    params.p    Number of lags
%    params.B    VAR coefficients
%    params.S    Covariance matrix
% =================================================================


F = randn(params.T+100,params.k);
% Normalize factors to be N(0,I) structural shocks
F = (grams(F')');   % First orthogonalize
F = F./std(F);      % Then standardize
e_t = F*params.L' + randn(params.T+100,params.M)*sqrt(params.R);
params.S = cov(e_t);
    
BETA  = params.B;
SIGMA = params.S;

%----------------------GENERATE--------------------------
% Set storage in memory for y
% First params.p rows are created randomly and are used as
% starting (initial) values
y =[0*rand(params.p,params.M) ; zeros(params.T+100-params.p,params.M)];

% Now generate Y from VAR (params.p,PHI,PSI)
for nn = params.p:params.T+100
    ylag = mlag2(y,params.p);
    y(nn,:) = ylag(nn,:)*params.B + e_t(nn,:);
end    
y = y(end-params.T+1:end,:);
