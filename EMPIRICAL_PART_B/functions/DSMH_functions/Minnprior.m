function [Mtildeinv,Sstar]  = Minnprior(yall,nlags,T,tstart,tend,n,r,lambda01,lambda11,lambda02,lambda12,lambda32)

%===============================================================
% calculate variance-covariance matrix of univariate AR residuals
e = zeros(T,n);
i = 0;
while i < n
   i = i+1;
   if nlags > 0
      ylags = yall(tstart-1:tend-1,i);
      ilags = 1;
      while ilags < nlags
         ilags = ilags + 1;
         ylags = [ylags yall(tstart-ilags:tend-ilags,i)];
      end
      ylags = [ylags ones(T,1)];
   else
      ylags = ones(T,1);
   end
   e(:,i) = yall(tstart:tend,i) - ylags*inv(ylags'*ylags)*ylags'*yall(tstart:tend,i);
end
Sstar = e'*e/T;

%============================================================
% Calculate inverse of Mtilde for standard Minn prior
v11 = 1:nlags;
v11 = v11'.^(-2*lambda11);
v21 = 1./diag(Sstar(1:r,1:r));
v31 = kron(v11,v21);
v31 = lambda01^2*v31;
v31 = 1./v31;

v12 = 1:nlags;
v12 = v12'.^(-2*lambda12);
v22 = 1./diag(Sstar(r+1:end,r+1:end));
v32 = kron(v12,v22);
v32 = lambda02^2*[v32; lambda32^2];
v32 = 1./v32;

Mtildeinv = diag([v31;v32]);