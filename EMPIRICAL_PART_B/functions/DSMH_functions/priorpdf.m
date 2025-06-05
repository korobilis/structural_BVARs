function [pdf_prior] = priorpdf(xbound,steps,cL,sigL,nuL,nL,signL)


% calculate prior densities to plot against posteriors
z1=-xbound:steps:xbound;           %grid for parameters to define points where density is evaluated
pdf_prior = zeros(nL,size(z1,2));

i = 0;
while i < nL
     i = i + 1;
     z2 = (z1 - cL(i))/sigL(i);
     if signL(i) == 1
 	     pdf_prior(i,:) = ((z1 > 0.0) .* (tpdf(z2,nuL(i))/sigL(i))) / (1 - tcdf(-cL(i)/sigL(i),nuL(i)));
     elseif signL(i) == -1
	     pdf_prior(i,:) = ((z1 < 0.0) .* (tpdf(z2,nuL(i))/sigL(i))) / tcdf(-cL(i)/sigL(i),nuL(i));
     elseif signL(i) == 0
	     pdf_prior(i,:) = tpdf(z2,nuL(i)) / sigL(i);
     end 
end