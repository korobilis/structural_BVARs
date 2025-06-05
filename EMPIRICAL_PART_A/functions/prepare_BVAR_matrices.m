function [Y,X,M,t,KK,K] = prepare_BVAR_matrices(Y,p)


% Number of observations and dimension of X and Y
t = size(Y,1); % t is the time-series observations of Y
M = size(Y,2); % M is the dimensionality of Y
KK = M*p+1;
K = KK*M;
	
% ===================================| VAR EQUATION |==============================
% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
t = t-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [ones(t,1) Ylag];
Y = Y(p+1:end,:);
