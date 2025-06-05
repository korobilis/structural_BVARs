
omegahat = (YY'*YY - YY'*XX*inv(XX'*XX)*XX'*YY)/T;
omegahat = eye(size(omegahat));

% prior for D
kappa = 2.0*ones(n,1);
% tau is determined from kappa and draw of A

% prior for B: partition in B1 and B2
%eta = 0.75*[eye(n) zeros(n,k-n)];
eta1 = 0.75*[eye(r) zeros(r,k-r)];
eta2 = zeros(n-r,k);
eta=[eta1;eta2];

% informativeness for B1
lambda01 = 0.5;   
lambda11 = 1.0;  
% informativeness for B2
lambda02 = 0.01;   
lambda12 = 1.0;  
lambda32 = 100.0;

    
%prior mode for elements in lambda
cL = [.8; .75; 0; 0; 0.1; .75; .75; 0; .7; .2; 0.6; -.1; -.05; ...
      -.5; 0; 0; 0; 0; .9; .9; 0; -.3; -.1; -.4; 0.2; 0; ... 
      -.3; -.15; 0.25; -.5; -.1; -.15; -.15; -.1; -.2; -.5; -.3; .25; -1];   % prior mode
sigL = 0.3*ones(nL,1);    % prior scale (determines how confident you are in your prior)
nuL = 3000*ones(nL,1);    % degrees of freedom: here set to approximate normal
%sign restrictions on impact effects where
%1=positive, -1=negative, 0=no restriction
signL = [1; 1; 0; 0; 1; 1; 1; 0; 1; 1; 1; -1; -1; ...
        -1; 0; 0; 0; 0; 1; 1; 0; -1; -1; -1; 1; 0; ...
        -1; -1; 1; -1; -1; -1; -1; -1; -1 ; -1; -1; 1; -1]; 