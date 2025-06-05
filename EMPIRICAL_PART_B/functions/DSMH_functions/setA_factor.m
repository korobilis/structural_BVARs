function A = setA_factor(vecL,n,r)
% input: vector of unique elements of loadings (vecL)
% output: full matrix A based on equation(2)
vecL1 = ones(n*r,1);
if n==r
   vecL(end+1,1) = 0;
   for jj = 1:r-1
       vecL1(jj*(n+1)-(n-1):jj*(n+1))=vecL(1+n*(jj-1):jj*n);
   end

   lam = reshape(vecL1,n,r);
   A = (lam'*lam)\lam';
        
elseif r<n
    for jj = 1:r-1
       vecL1(jj*(n+1)-(n-1):jj*(n+1))=vecL(1+n*(jj-1):jj*n);
    end
    vecL1(r*(n+1)-(n-1):end)=vecL(1+n*(r-1):end);
    lam = reshape(vecL1,n,r);
    H2 = [zeros(n-r,r) eye(n-r)]*(eye(n)-lam*((lam'*lam)\lam'));
    H2_tilde = zeros(size(H2));
    H2_tilde(1,:) = H2(1,:);
    
    for jx = 2:(n-r)
        y = H2(jx,:)';
        x = H2(1:jx-1,:)';
        beta = x'*x\x'*y;
        H2_tilde(jx,:) = (y - x*beta)';
    end
    
    A = [(lam'*lam)\lam';H2_tilde];
end