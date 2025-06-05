function [ytilde,yxtilde,xtildei,kappastar] = posteriorquantities(YY,XX,Mtildeinv,eta,k,n,kappa,T,prior)

%============================================================
% calculate parameters that characterize posteriors

kappastar = kappa + (T/2); 
if prior==0
    ytilde = YY'*YY + eta*Mtildeinv*eta';
    yxtilde = YY'*XX + eta*Mtildeinv;
else
    ytilde = eta*Mtildeinv*eta';
    yxtilde = eta*Mtildeinv;
end
    
xtildei = zeros(k,k,n);
i = 0;
while i < n
    i = i+1;
        if prior == 0
            xtildei(:,:,i) = inv(XX'*XX + Mtildeinv);
        else 
            xtildei(:,:,i) = inv(Mtildeinv);
        end
end