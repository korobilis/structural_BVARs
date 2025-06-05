n1 = 3;
n2 = 5;

% anames = {' \pi^{d}';' ffr^{d}'; ' spread^{d}'; ' comm^{d}'; ' wage^{d}'; ...
%           ' gap^{s}';' ffr^{s}'; ' spread^{s}'; ' comm^{s}'; ' wage^{s}'; ...
%           ' gap^{mp}';' \pi^{mp}'; ' spread^{mp}'; ' comm^{mp}'; ' wage^{mp}'};

figure(2)
nbin=50;
i = 0;
while i < nL
     i = i+1;
     %subplot(n1,n2,i)
     [ag,bg]=hist(allgDraws(:,i),nbin);
     delta=bg(1,2)-bg(1,1);
     post_a = ag./((10000)*delta);
     bar(bg,post_a), hold on, plot(z1,pdf_prior(i,:),'r','linewidth',2); box on
     %axis([-xbound xbound 0 ybound])
     %title_a = strcat(anames(i));
     %title(title_a,'fontsize',12)
     %pause
     
end
