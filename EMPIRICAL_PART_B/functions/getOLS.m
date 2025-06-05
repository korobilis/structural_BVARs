function [beta_OLS,beta_OLS_vec,SIGMA_OLS] = getOLS(y_t,x,M,p)

% Get OLS estimates
beta_OLS=((x'*x)\(x'*y_t));
bb1=beta_OLS(1,:);
bb2 = beta_OLS(2:end,:);
bb3 =[];
for i = 1:p
    temp = bb2((i-1)*M+1:i*M,:);
    bb3 = [bb3 ; temp(:)];
end
beta_OLS_vec = [bb1';bb3];
SIGMA_OLS = (y_t-x*beta_OLS)'*(y_t-x*beta_OLS)./(size(y_t,1)-p-1);