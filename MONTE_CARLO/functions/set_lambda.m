function L=set_lambda(vecL,n,r)

vecL1 = ones(n*r,1);
for jj = 1:r-1
    vecL1(jj*(n+1)-(n-1):jj*(n+1))=vecL(1+n*(jj-1):jj*n);
end
vecL1(r*(n+1)-(n-1):end)=vecL(1+n*(r-1):end);

L = reshape(vecL1,n,r);