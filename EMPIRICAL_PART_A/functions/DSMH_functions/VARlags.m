function [XX,T,k,n,nL] = VARlags(YY,yall,nlags,r,tstart,tend)

% ==================================================
% construct matrix of lags 
%      XX = (T x k) matrix of observations on k different regressors
%    (this section should work as is for all applications)

T = size(YY,1);     % T is sample size
n = size(YY,2);     % n is number of equations/variables
nL = (n-1)*r;       % number of elements in L for which priors are formed

if nlags > 0
    XX = yall(tstart-1:tend-1,:);
    ilags = 1;
    while ilags < nlags
        ilags = ilags + 1;
        XX = [XX yall(tstart-ilags:tend-ilags,:)];
    end
    XX = [XX ones(T,1)];
else
    XX = ones(T,1);     % always include constant in XX 
end

k = size(XX,2);         % k is number of regressors