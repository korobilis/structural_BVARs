function [Beta] = noninformative(y,X,sigma_sq,type,ieq)


[n,p]=size(X);

%% paramters %%
Beta=zeros(p,1); 

%% matrices %%
I_n=eye(n); 
l=ones(n,1);
if type==2
   Q_star=X'*X;
end

lambda_star = 1*ones(p,1);
lambda_star(1) = 100;
lambda_star(1+ieq) = 1;

%% update beta %%
switch type

    case 1 % new method
        U=bsxfun(@times,(lambda_star.^2),X');
        %% step 1 %%
        u=normrnd(0,lambda_star);
        v=X*u + randn(n,1);
        %% step 2 %%
        v_star=(X*U+I_n)\((y/sqrt(sigma_sq))-v);
        Beta=sqrt(sigma_sq)*(u+U*v_star);
        
    case 2 % Rue  
        L=chol((1/sigma_sq)*(Q_star+diag(1./lambda_star.^2)),'lower');
        v=L\((y'*X)./sigma_sq)';
        mu=L'\v;
        u=L'\randn(p,1);
        Beta=mu+u;  
end

