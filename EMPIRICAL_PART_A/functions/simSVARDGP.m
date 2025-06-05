function [y, BETA, SIGMA, C, R]   =  simSVARDGP(DGP,params)

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

if nargin == 1
    params.T = 200;
    params.M = 12;
    params.k = 3;
    params.p = 1;
  
    % Note: we allow AR(1) value to be random in the range (0.4, 0.6) and then set for a certain VAR equation to have
    % persistence: ar(1)/(1^2) + ar(1)/(2^2) + ar(1)/(3^2) + ar(1)/(4^2).... in the spirit of the Minnesota prior. For
    % ar(1) = 0.6 we get a total persistence of 0.8542 for a 4-lag specification (while asymptotically the persistence
    % converges to ~0.987). Notice that I allow the persistence to be random for each VAR equation
    % Regarding correlations, I work with the matrix C where I randomly set some lower triangular elements to zero and
    % some to a value between (0,1). If we don't use SSVS in the comparison for this DGP, then I think it would be ok to
    % define things in terms of C.    
    low_lim = 0.7;   high_lim = 0.8;
    lam     = (2./(params.M-1)); % this makes the prob of non-zero elements in the other lags a function of the VAR size; 
                                 % N=3 => lam=0.5; N=7 => lam=0.17; N=20 => lam=0.05; N=40 => lam=0.025; ....  
    sig2_sl = 0.05;              % variance of non-zero coeffs in the other lags
    
    stable_var = 0;
    while stable_var == 0
        persistence_level = (high_lim - low_lim)*rand(params.M,1) + low_lim;
                
        % Fill in own lags                
        PHI = ones(1,params.M);   % Constant
        for i = 1:params.p
            PHI = [PHI; (diag(persistence_level)./(i^2))]; 
        end
        
        % Fill in other lags
        PHI_long = PHI(:);
        indx_other = PHI_long==0;
        PHI_long(indx_other) = 1*double(rand(sum(indx_other),1)) .* (sqrt(sig2_sl)*randn(sum(indx_other),1));
%         PHI_long(indx_other) = (sqrt(sig2_sl)*randn(sum(indx_other),1));
        params.B = reshape(PHI_long,params.M*params.p+1,params.M);
        
        % Check whether VAR is stable, writing VAR in companion form
        By = [params.B(2:end,:)'; eye(params.M*(params.p-1)) , zeros(params.M*(params.p-1),params.M)];
        if max(abs(eig(By)))>=1
            disp('Non stationary VAR - redrawing now...')
            stable_var = 0;
        else
            stable_var = 1;
        end
        
        % Correlation and Covariance matrices (C has lower-dimensional structure)
        params.cL = [0.8; 0.75; 0.15; 0.25; 0.05; -0.5; 0.25; -0.9; 0.6; 0.15; -0.3; 0.8; 0.25; -0.5; -0.1; 0.2; ...
            0.3; 0.75; 0.5; 0.4; 0.1; -0.2; 0.9; 0.5; 0.2; 0.4; -0.6; -0.2; 0.7; -0.25; -0.4; 0.3; -0.5];
        C = set_lambda(params.cL,params.M,params.k);
        R = 0.01*diag(rand(params.M,1));
        if DGP == 1
            S = eye(params.k);%diag(rand(params.k,1));            
            params.S = C*S*C' + R;
        elseif DGP == 2
            F = randn(params.T+1000,params.k);
            % Normalize factors to be N(0,I) structural shocks
            F = (grams(F')');   % First orthogonalize
            F = F./std(F);      % Then standardize
            e_t = F*C' + randn(params.T+1000,params.M)*sqrt(R);
            params.S = cov(e_t);
        end
    end
end
    
BETA  = params.B;
SIGMA = params.S;

%----------------------GENERATE--------------------------
% Set storage in memory for y
% First params.p rows are created randomly and are used as
% starting (initial) values
y =[0*rand(params.p,params.M) ; zeros(params.T+1000-params.p,params.M)];

% Now generate Y from VAR (params.p,PHI,PSI)
if DGP == 1
    for nn = params.p:params.T+1000
        u = chol(params.S)'*randn(params.M,1);
        ylag = mlag2(y,params.p);
        y(nn,:) = [1 ylag(nn,:)]*params.B + u';
    end    
elseif DGP == 2
    for nn = params.p:params.T+1000
        ylag = mlag2(y,params.p);
        y(nn,:) = [1 ylag(nn,:)]*params.B + e_t(nn,:);
    end    
end
y = y(end-params.T+1:end,:);
