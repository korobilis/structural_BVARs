function p = likeli_val(theta,kappa,T,omegahat,Sstar,...
    ytilde,xtildei,yxtilde,r)
% ===================================================
% By Christiane Baumeister and James D. Hamilton (Nov 2014)
% inputs
   % theta = value of the (nA x 1) vector of unknown elements in matrix A
   % A_params = (nA x 4) matrix of parameters characterizing priors for elements of theta
   % kappa = (n x 1) vector of parameters characterizing priors for D
   % T = sample size
   % kexample governs how theta maps into matrix A
   % omegahat = (n x n) matrix of OLS residual covariances
   % Sstar = (n x n) matrix of univariate residual covariances 
   % Ri, Vi = params for long-run priors (set to zero if not used) 
% outputs
  % p = negative of value of log posterior distribution evaluated at the point theta
% ===================================================
% evaluate log of prior for A at theta
   n = size(omegahat,1);
      
% evaluate log of other elements of posterior
   A = setA_factor(theta,n,r);
   Q = A*omegahat*A';
   tau = kappa.*diag(A*Sstar*A');
   kappastar = kappa + (T/2);
   zetastar = zeros(n,1);
   i = 0;
   while i < n
       i = i+1;
       ytildei = A(i,:)*ytilde*A(i,:)';
       yxtildei = A(i,:)*yxtilde;
       zetastar(i) = ytildei - yxtildei*xtildei(:,:,i)*yxtildei';
   end
   
   
  p = (T/2)*log(det(Q)) - kappastar'*log((2*tau/T) + zetastar/T) + kappa'*log(tau);
  
   
