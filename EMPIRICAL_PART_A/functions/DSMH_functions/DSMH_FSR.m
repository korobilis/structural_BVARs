function [Bdraw,invdDraw,L_save,psi_invA0] = DSMH_FSR(yall,nlags,sg,hmax)

% Bayesian structural VAR with factor model contemporaneous decomposition of the
% covariance matrix and sign restrictions
%
% This code uses the DSMH algorithm in Baumeister and Korobilis (2020)
%-------------------------------------------------------------------------------
% Written by Christiane Baumeister & Dimitris Korobilis
% University of Notre Dame & University of Glasgow
% This version: 10 February 2020
%-------------------------------------------------------------------------------


%-------------------------------PRELIMINARIES--------------------------------------
% User-defined settings for the DSMH algorithm
nH=30;           %number of tempering stages
nM=20;           %number of striations
nN=500;          %length of Markov chains in each stage
nG=20;           %number of parallel Markov chains in each stage
nK=500;          %length of the training chain to determine ci
Tau=5;           %thinning parameter (only every Tau-th draw of the Markov chain is stored)
alpha0=0.2;      %lower bound for acceptance rate
alpha1=0.3;      %upper bound for acceptance rate
lambda_1=0.005;  %starting value of the tempering schedule

prior=0;   % set prior=1 to compute the IRFs implied by the prior
           % set prior=0 to compute the posterior IRFs

%% ===================================| SVAR EQUATION |==============================
[YY,XX,T,k,n] = prepare_SVAR_matrices(yall,nlags);

%% ==================| 2. PRIORS:
% Stuff for impulse response analysis: 
r = size(sg,2);        % Number of factors/structural shocks
sg(logical(eye(size(sg)))) = [];
signL = sg';
%signL = [1;0;0;0;0;0;-1;-1;1;1;1;1;1;-1;0;0;0;0;-1;-1;1;-1;-1;1;1;1;-1;-1;-1;0;-1;-1;-1;1;-1;-1;-1;-1;-1];
nL = (n-1)*r; % Number of unrestricted elements in L (we have r diagonal coefficients fixed to 1)

% set values that characterize the priors and how they will be plotted
omegahat = (YY'*YY - YY'*XX*inv(XX'*XX)*XX'*YY)/T;
omegahat = eye(size(omegahat));

% prior for D:
kappa = 2.0*ones(n,1); % tau is determined from kappa and draw of A

% prior for B: 
% partition in B1 and B2
eta1 = 1*[eye(r) zeros(r,k-r)];    
eta2 = zeros(n-r,k);      eta=[eta1;eta2];
% informativeness for B1
lambda01 = 0.5;  lambda11 = 1.0;  
% informativeness for B2
lambda02 = 0.01;    lambda12 = 1.0;   
lambda32 = 100.0;    
 % Get here Minnesota prior moments for SVAR
[Mtildeinv,Sstar]  = Minnprior(yall,nlags,T,nlags+1,T+nlags,n,r,lambda01,lambda11,lambda02,lambda12,lambda32);

% prior for L: 
% prior mode for elements in lambda
cL = [1;0.5;0;-0.5;-0.7;0;-0.7;-0.6;1;0.8;1.2;0.7;0.6;-1.5; 0;0; 0; 0;-0.5;-0.5;0.5;...
    -1;-1;1;0.5;0.5;-0.50;-1.5;-1;-3;-2;-0.44;-2;0.2;-0.5;-0.5;-2;-2.5;-1.5];
if n*r < length(cL)+r
    Ltemp = set_lambda(cL,14,3); Ltemp = Ltemp(1:n,1:r); Ltemp(logical(eye(size(Ltemp))))=[]; cL = Ltemp';
elseif n*r > length(cL)+r
    Ltemp = zeros(n*r-r,1); Ltemp(1:39) = cL; cL = Ltemp;
end
sigL = 1*ones(nL,1);      % prior scale (determines how confident you are in your prior)
nuL = 3000*ones(nL,1);    % degrees of freedom: here set to approximate normal

% Evaluate prior/posterior
xbound = 2;    ybound = 2;   steps = 0.01;
[pdf_prior] = priorpdf(xbound,steps,cL,sigL,nuL,nL,signL);
[ytilde,yxtilde,xtildei,kappastar] = posteriorquantities(YY,XX,Mtildeinv,eta,k,n,kappa,T,prior);

%% DSMH setup %%
f_anon = @(theta)likeli_val(theta,kappa,T,omegahat,Sstar,ytilde,xtildei,yxtilde,r);
c=size(cL,1);   
%tempering schedule (geometric progression)
g=(1:1:nH)';
lambda=zeros(nH,1);
lambda(1,1)=1/exp(50);
for jx=2:nH
    lambda(jx,1)=lambda_1.^((nH-g(jx,1))/(nH-1));
end
p=1/(20*Tau);   % probability that the proposal comes from the distribution at the previous stage
nStr=nG*nN/nM;  % striation size

%create a population of draws from joint prior distribution of model parameters and evaluate 
MCsample=[];
MCeval=[];
parfor jx=1:nN*nG
    [a,b]=get_prior_draw(c,cL,sigL,nuL,signL);
    MCsample=[MCsample;a];
    MCeval=[MCeval;b];
end

[valp,indx]=sort(MCeval);   %get index of ordered values
theta=MCsample(indx,:);     %order the draws wrt target fct to build striations

weights=ones(nN*nG,1)./(nN*nG);          %importance weights
ci=1;                                    %tuning parameter for MH step
wtheta=weights.*theta;                   
Sigi=(theta'*wtheta)-(wtheta'*wtheta);   %variance matrix for MH step

allgDraws=zeros(nN*nG,c);
allgDrawVals=zeros(nN*nG,1);

%% here starts the big s-th stage loop
tic;
for s=2:nH
    disp(['s = ' num2str(s)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sub-routine to determine ci at stage s (see Appendix B)
    StartPars=datasample(theta,nG,'Weights',weights,'Replace',false);
    lam=lambda(s,1);
    nacc=[];
    alpha_c=0;
    
    while alpha_c == 0
        parfor jj=1:nG
            old = StartPars(jj,:);          
            ptarget_old = logP_factor(old',cL,sigL,nuL,signL) + lam*f_anon(old');
            naccept = 0;                             
            for ii=1:nK
                proposal = old + (chol(ci*Sigi)'*randn(c,1))';           
                ptarget_new = logP_factor(proposal',cL,sigL,nuL,signL) + lam*f_anon(proposal');
            
                accept=min([exp(ptarget_new - ptarget_old);1]);
                if  rand <= accept
                    old = proposal;                 %we retain the new draw
                    ptarget_old = ptarget_new;
                    naccept = naccept + 1;          %count the number of acceptances
                end                
            end            
            nacc=[nacc;naccept];
        end
        
        alpha_c = sum(nacc)/(nK*nG);        
        if alpha_c < alpha0 || alpha_c > alpha1
            ci = ci*CiMultiple(alpha_c,alpha0,alpha1);
            alpha_c = 0;
            nacc = [];
        end        
    end
    disp(alpha_c)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run striated MH sampler at stage s
    sindx=datasample((1:1:nN*nG)',nG,'Weights',weights,'Replace',false);
    StartPas=theta(sindx,:);
    lamp=lambda(s-1,1);  %previous lambda
    
    gDraws=[];
    gDrawVals=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5   
    parfor jx=1:nG
        old = StartPas(jx,:);
        ptarget_old = logP_factor(old',cL,sigL,nuL,signL) + lam*f_anon(old');
           
        jj=1;        
        while jj<=(nN*Tau)
            if rand < p   %the draw comes from the striation                   
                stria = ((ceil(sindx(jx,1)/nStr)-1)*nStr+1:1:ceil(sindx(jx,1)/nStr)*nStr)';
                Strindx = datasample(stria,1,'Replace',false);
                proposal = theta(Strindx,:);
                ptarget_new = logP_factor(proposal',cL,sigL,nuL,signL) + lam*f_anon(proposal');
            
                fn_s_old = logP_factor(old',cL,sigL,nuL,signL) + lamp*f_anon(old');            
                fn_s_proposal = valp(Strindx,1);                   
                accept = min([exp(ptarget_new - ptarget_old + fn_s_old - fn_s_proposal);1]);                
                if rand <= accept
                    old = proposal;                  %we retain the new draw        	       
                    ptarget_old = ptarget_new;
                end
                jj=jj+1;
                cc=1;                              
            else   %the draw comes from the regular MH move                   
                proposal = old + (chol(ci*Sigi)'*randn(c,1))';                
                if min(sign(proposal').*signL) >=0                       
                   ptarget_new = logP_factor(proposal',cL,sigL,nuL,signL) + lam*f_anon(proposal');                   
                   accept = min([exp(ptarget_new - ptarget_old);1]);
                   if rand <= accept
                       old = proposal;                %we retain the new draw
                       ptarget_old = ptarget_new;
                   end
                   jj=jj+1;
                   cc=1;                                     
                else
                    cc=0;
                end                
            end
            if cc==1
                gDraws = [gDraws;old];
                gDrawVals = [gDrawVals;exp(ptarget_old)];
            end            
        end           
    end    
    counter=0;
    for jjj=5:5:(nN*Tau*nG)
        counter                 = counter + 1;
        allgDraws(counter,:)    = gDraws(jjj,:);       
        allgDrawVals(counter,1) = gDrawVals(jjj,1);          
    end
    
    weight_unnorm = zeros(nN*nG,1);
    parfor jy=1:nN*nG
        draw = allgDraws(jy,:);
        denom = exp(logP_factor(draw',cL,sigL,nuL,signL) + lamp*f_anon(draw'));
        weight_unnorm(jy,1) = allgDrawVals(jy,1)/(denom  + 1e-9);
    end
    weights = (weight_unnorm + 1e-9)./sum(weight_unnorm + 1e-9);
    [valp,indx]=sort(allgDrawVals);
    theta=allgDraws(indx,:);
    
    wtheta=weights.*theta;                   
    Sigi=(theta'*wtheta)-(wtheta'*wtheta);   %variance matrix for MH step    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;

L_save = [];
for i = 1:r 
    L_save = [L_save, [ones(counter,1), allgDraws(:,(i-1)*14+1:min(i*14,r*(n-1))) ]];
end
L_save = reshape(L_save,counter,n,r);

%figure2
nuse   = size(allgDraws,1);
a_post = allgDraws';

for jx=1:size(a_post,2)
   A = setA_factor(a_post(:,jx),n,r);
   Q = A*omegahat*A';
   tau = kappa.*diag(A*Sstar*A');
   i = 0;
   while i < n
       i = i+1;
       ytildei = A(i,:)*ytilde*A(i,:)';
       yxtildei = A(i,:)*yxtilde;
       zeta_post(i,jx) = ytildei - yxtildei*xtildei(:,:,i)*yxtildei';
       mstar_post(i,:,jx) = yxtildei*xtildei(:,:,i);
   end
end

%%% This program generates draws for D and B as well as IRFs
kcum = 0;   % 0: not cumulated; 1: cumulated for IRFs and var decomps

cholstar_i = zeros(k,k,n);  % cholstar_i(:,:,j) is Cholesky factor of Mstar(j)
j = 0;
while j < n
    j = j+1;
    cholstar_i(:,:,j) = chol(xtildei(:,:,j))';
end

% create blank matrices in which posterior draws will be filed in
invdDraw = zeros(n,nuse);   % draws for diag(inv(D)) will go in columns
Bdraw = zeros(n,k,nuse);    % draws for B go in as matrices

% ===========================================================
% Generate draws for D and B
disp('Generating draws for D and B')
isim = 0;
while isim < nuse
    isim = isim + 1;
    A = setA_factor(a_post(:,isim),n,r);   
    % generate a draw for diag(inv(D)), put in invdDraw
    taustar = kappa.*diag(A*Sstar*A') + zeta_post(:,isim)/2;
    invdDraw(:,isim) = gamrnd(kappastar,1./taustar);
    
    % generate a draw for (n x k) matrix B, put in BDraw
    j = 0;
    while j < n
        j = j+1;
        Bdraw(j,:,isim) = mstar_post(j,:,isim) + ...
            (1/sqrt(invdDraw(j,isim)))*(cholstar_i(:,:,j)*randn(k,1))';
    end
    
    if (isim/10000) == floor(isim/10000)
        disp(isim)
    end
    
end

% ====================================================================
% calculate impulse-response function for each draw
hsim = 0;    % hsim is horizon being calculated
if nlags == 0
    hmax = 1;
end

psi = zeros(n,n,hmax+nlags,nuse);  
    % psi(:,:,hsim,isim) is n x n matrix of nonorthogonalized IRF
    % first nlags -1 elements are defined to be zero to use same recursion
       % for horizon hsim and draw isim
psi_invA = zeros(n,n,hmax+nlags,nuse);
    % psi_invA(:,:,hsim,isim) is n x n matrix of structural IRF
psi_invA_cum = zeros(n,n,hmax+nlags,nuse);
    % this collects cumulated IRFs if desired
phicheck = zeros(n,k);
countA = zeros(n,n,nlags+hmax); % countA(:,:,hsim) is fraction of positive draws
isim = 0;
disp('Calculating impulse-response functions')
while isim < nuse
    isim = isim+1;
    A = setA_factor(a_post(:,isim),n,r);
    invA = inv(A);
    psi(:,:,nlags,isim) = eye(n);
    psi_invA(:,:,nlags,isim) = invA;
    psi_invA_cum(:,:,nlags,isim) = invA;
    countA(:,:,nlags) = countA(:,:,nlags) ...
        + (psi_invA(:,:,nlags,isim) > zeros(n,n));
    phi = invA*Bdraw(:,:,isim);
    phi_sim(:,:,isim) = phi;
    phicheck = phicheck + phi;
    hsim = 1;
   
    while hsim < hmax
        ilags = 0;
        while ilags < nlags
            ilags = ilags+1;
            psi(:,:,nlags+hsim,isim) = psi(:,:,nlags+hsim,isim) ...
               + phi(:,(ilags-1)*n+1:ilags*n)*psi(:,:,nlags+hsim-ilags,isim);
        end
        psi_invA(:,:,nlags+hsim,isim)=psi(:,:,nlags+hsim,isim)*invA;
        psi_invA_cum(:,:,nlags+hsim,isim)=psi_invA_cum(:,:,nlags+hsim-1,isim)+...
                psi(:,:,nlags+hsim,isim)*invA;
        countA(:,:,hsim+nlags) = countA(:,:,hsim+nlags) ...
            + (psi_invA(:,:,nlags+hsim,isim) > zeros(n,n));
        hsim = hsim+1;
    end
    
    if (isim/10000) == floor(isim/10000)
        disp(isim)
    end
end


% ===================================================================
% plot graphs of impulse response function
index1 = round(0.025*nuse);
index2 = round((1 - 0.025)*nuse);
index3 = round(0.16*nuse);
index4 = round((1 - 0.16)*nuse);
HO=(0:1:20)';
if kcum == 0
    psi_invA0 = psi_invA;
elseif kcum == 1
    psi_invA0 = psi_invA_cum;
end

