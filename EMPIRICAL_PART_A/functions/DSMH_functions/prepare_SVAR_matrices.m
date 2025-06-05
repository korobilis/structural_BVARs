function [Y,X,t,k,n] = prepare_SVAR_matrices(Y,p)


% Number of observations and dimension of X and Y
t = size(Y,1); % t is the time-series observations of Y
n = size(Y,2); % M is the dimensionality of Y
k = n*p+1;

% ===================================| VAR EQUATION |==============================
% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
t = t-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [Ylag ones(t,1)];
Y = Y(p+1:end,:);


