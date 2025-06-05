function [beta_save,SIGMA_save,L_save,F_save,irf_sign,DIC,DIC2] = Gibbs_FSR_TIMES(Y,p,sg,nhor,ngibbs,nburn,nthin)

[y,x,M,T,KK,K] = prepare_BVAR_matrices(Y,p);

[beta_OLS,~,~] = getOLS(y,x,M,p);

sg0 = sg; 
sg0(isnan(sg0)) = 0;
nor = numel(find(sg~=0));  % number of restrictions 
N   = size(sg,2); % number of structural shocks

% Priors on model parameters
if size(x,1)>size(x,2)
    est_alg = 2;
else
    est_alg = 1;
end

psi = cell(M,1); tau = zeros(M,1); sigma_sq = zeros(M,1);
for ieq = 1:M
    psi{ieq,1} = ones(KK,1);
    tau(ieq,1) = 1;
    sigma_sq(ieq,1) = 1;
end

lam_prmean  = zeros(M,N);
lam_prvar   = 4*ones(M,N);
lam_iprvar  = 1./lam_prvar;
index_vec   = 1:N;

% Extract factors
yhat = y - x*beta_OLS;
[F,L] = extract(zscore(yhat),N);

%========|STORAGE MATRICES:
R = eye(M);
iR = ones(M,1);
beta_mat = zeros(KK,M);

beta_save  = zeros(ngibbs/nthin,KK,M);
SIGMA_save = zeros(ngibbs/nthin,M,M);
L_save     = zeros(ngibbs/nthin,M,N);
F_save     = zeros(ngibbs/nthin,T,N);
R_save     = zeros(ngibbs/nthin,M);

irf_sign   = zeros(ngibbs/nthin,M,N,nhor);

D_theta    = zeros(ngibbs/nthin,1);
D_theta2   = zeros(ngibbs/nthin,1);
%======================= BEGIN MCMC ESTIMATION =======================
tic;
counter = 0;
for iter = 1:ngibbs + nburn
    if mod(iter,1000)==0   
        disp(['This is iteration ' num2str(iter)])         
        disp([num2str(100*(iter/(ngibbs+nburn))) '% completed'])
        toc;
    end
    
    y_til = y - F*L';
    % STEP 1: update BETA
    for ieq = 1:M
       [beta_i,psi{ieq,1},tau(ieq)] = horseshoe(y_til(:,ieq),x,psi{ieq,1},tau(ieq),R(ieq,ieq),est_alg,ieq);
       beta_mat(:,ieq) = beta_i;         
    end
    B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
    
%     % check BETA for stationarity
%     rej_count = 0;
%     while max(abs(eig(B))) > 0.999
%         rej_count = rej_count + 1;
%         if rej_count > 100; y_til  = y; end
%         for ieq = 1:M
%             [beta_i,psi{ieq,1},tau(ieq)] = horseshoe(y_til(:,ieq),x,psi{ieq,1},tau(ieq),R(ieq,ieq),est_alg,ieq);
%             beta_mat(:,ieq) = beta_i;         
%         end
%         B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
%     end
    
    % STEP 2: update SIGMA 
    yhat = y - x*beta_mat;  % yhat is the VAR residual that follows the factor model          

    % Sample shocks
    for t = 1:T
        F_var  = inv(1*eye(N) + L'*diag(iR)*L);
        F_mean = F_var*L'*diag(iR)*yhat(t,:)';
        F(t,:) = (F_mean + chol(F_var)'*randn(N,1))';
    end
	F = F - mean(F); % demean factor sample
    
    % I am doing now a very inefficient element-by-element update of the posterior of lambda,
    % so I can make a straightfoward use of the univariate truncated normal generator 
    for i = 1:M        
        %L_var_post  = inv( lam_iprvar + (F'*F).*iR(i) );
        %Lpostmean   = L_var_post*((F'*yhat(:,i)).*iR(i));
        Lvar = diag(lam_iprvar(i,:)) + (F'*F).*iR(i); % inverse of posterior covariance matrix of L(i,:)
        L_bar(i,:) = Lvar\(diag(lam_iprvar(i,:))*lam_prmean(i,:)' + (F'*yhat(:,i)).*iR(i));
        for j = 1:N
            whole_vec = Lvar(j,:)'.*(L(i,:) - L_bar(i,:));
            Lpostvar = 1./Lvar(j,j); ss = sqrt(Lpostvar);
            Lpostmean = L_bar(i,j) - sum(whole_vec(find(index_vec~=j)))*Lpostvar;
            if sg(i,j)  == 1
                lower  = (0 - Lpostmean)/ss;   
                upper  = (10000 - Lpostmean)/ss;
                L(i,j) = trandn(lower,upper);
                L(i,j) = Lpostmean + ss*L(i,j);
            elseif sg(i,j)  == -1
                lower  = (-10000 - Lpostmean)/ss;
                upper  = (0 - Lpostmean)/ss;
                L(i,j) = trandn(lower,upper);
                L(i,j) = Lpostmean + ss*L(i,j);
            elseif sg(i,j) == 0
                L(i,j)  = 0;
            elseif isnan(sg(i,j))
                L(i,j)  = Lpostmean + ss*randn;
            end
        end
    end
%     % Normalize diagonal of loadings to be ones
%     L = L/diag(diag(L));
    
    % draw R
    sse2    = sum((yhat - F*L').^2);
    R_1     = (1 + T);    
    R_2     = (0.1 + sse2);
    iR      = gamrnd(R_1./2,2./R_2);
    R       = diag(1./iR);
    
    SIGMA  = L*L' + R;

    if iter>nburn && mod(iter,nthin) == 0
        counter = counter + 1;
        beta_save(counter,:,:)  = beta_mat; 
        SIGMA_save(counter,:,:) = SIGMA;
        L_save(counter,:,:)     = L;
        F_save(counter,:,:)     = F;
        R_save(counter,:)       = diag(R);
        
        % ------------------| IMPULSE RESPONSE ANALYSIS |------------------
        % 1) Obtain propagation coefficients in companion form       
        B = [beta_mat(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
        shock = L;
                       
        % 2) Obtain IRFs
        IRF = zeros(M,N,nhor);
        for ihor = 1:nhor
            Blarge = B^(ihor-1);
            IRF(:,:,ihor) = Blarge(1:M,1:M)*shock;
        end
        
        % save draws
        irf_sign(counter,:,:,:) = IRF;
        
        % Calculate DIC        
        c = -T*M/2*log(2*pi);
        u = y - x*beta_mat - F*L'; u = reshape(u',T*M,1);
        D_theta(counter,:)  = c - T/2*log(det(R)) - .5*u'*kron(speye(T),diag(iR))*u;
        D_theta2(counter,:) = -T*8/2*log(2*pi) - T/2*log(det(R(1:8,1:8))) - .5*u(1:T*8)'*kron(speye(T),diag(iR(1:8)))*u(1:T*8);
    end
end
%======================== END MCMC ESTIMATION ========================
BB = squeeze(median(beta_save,1));
LL = squeeze(median(L_save,1));
RR = squeeze(median(R_save,1));  iRR = diag(1./RR);
FF = squeeze(median(F_save,1));

u = y - x*BB - FF*LL'; u = reshape(u',T*M,1);
D_theta_bar  = c - T/2*log(det(diag(RR))) - .5*u'*kron(speye(T),iRR)*u;
D_theta_bar2 = -T*8/2*log(2*pi) - T/2*log(det(diag(RR(1:8)))) - .5*u(1:T*8)'*kron(speye(T),iRR(1:8,1:8))*u(1:T*8);

DIC  = -4*mean(D_theta) + 2*D_theta_bar;
DIC2 = -4*mean(D_theta2) + 2*D_theta_bar2;

