%%% This program generates draws for D and B as well as IRFs
%%% By Christiane Baumeister and James D. Hamilton (Nov 2014)

% To be run after running metrop_hast, which should have set variables including:
%  a_post = (nA x nsim) matrix; col isim is a simulated draw for vecstar(A)
%  zeta_post = (n x nsim) matrix; col isim is simulated draw for zetastar
%  mstar_post = (n x k x nsim);  mstar_post(i,:,isim) is (1 x k) vector of draws
%                                for coeffs of ith structural equation
% xtildei = (k x k x n);   xtildei(:,:,i) is (k x k) matrix of Mstar(i)
%  kexample to control call to priors and definition of A
%  nlags = number of lags in VAR
% before calling set kcum = 1 to accumulate IRF, = 0 for no accumulate

%=============================================================
% set initial parameters
if (ndraws-nburn)>500000
    nuse=500000;      %fewer draws because of memory constraints
else nuse=size(a_post,2);
end

% nsim = size(a_post,2);
% nuse = nsim;   % nuse is number of sims to be used, set = nsim for final

% set parameters that will be held constant for all draws
 
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
'Generating draws for D and B'
isim = 0;
while isim < nuse
    isim = isim + 1;
    A = setA(a_post(:,isim),n);   
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
        isim
    end
    
end

if n==6 && nuse>=1e6
    clear mstar_post zeta_post   %to free up space
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
'Calculating impulse-response functions'
while isim < nuse
    isim = isim+1;
    A = setA(a_post(:,isim),n);
    invA = inv(A);
    psi(:,:,nlags,isim) = eye(n);
    psi_invA(:,:,nlags,isim) = invA;
    psi_invA_cum(:,:,nlags,isim) = invA;
    countA(:,:,nlags) = countA(:,:,nlags) ...
        + (psi_invA(:,:,nlags,isim) > zeros(n,n));
    phi = invA*Bdraw(:,:,isim);
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
        isim
    end
end

ntable = 0;
ntablemax = 3;
'Fraction of responses greater than zero at horizon:'
while ntable < ntablemax
    ntable
    countA(:,:,ntable+nlags)/nuse
    post_probs(:,:,ntable+1)=countA(:,:,ntable+nlags)/nuse;
    ntable = ntable + 1;
end

if n==3
    save posterior_probs_table post_probs 
end
%'check: avg phi was'
%phicheck/nuse

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

if n==3
   load prior_IRFs
   figure(4)
   for ii = 1:n
       for jj = 1:n
           subplot(n,n,jj+(ii-1)*n)
           plotPanel_posterior;
       end
   end

elseif n==6
    figure(9)
    for ii = 1:3
        for jj = 1:3
            subplot(3,3,jj+(ii-1)*3)
            plotPanel_posterior6;
        end
    end
    figure(10)
    for ii = 1:n
        for jj = 3:3
            subplot(3,2,ii)
            plotPanel_posterior6;
        end
    end
     
end


