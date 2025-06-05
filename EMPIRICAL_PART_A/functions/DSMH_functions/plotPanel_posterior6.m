title_a = strcat('Response of',varnames(ii),' to', shocknames(jj));
x1 = (squeeze(psi_invA0(ii+15,jj,nlags:nlags+hmax-1,:)));  %%% add +3,+6 etc to ii to plot the next 3x3 plot of IRFS
x1 = sort(x1,2);
temp1=[(median(x1,2)) x1(:,index3) x1(:,index4)];   % posterior median and 68% credibility set
temp1=temp1(1:21,:);

plotx1(temp1,HO); box on; plot(HO,zeros(21,1),'k:','linewidth',1)
% hold on, plot(HO,x1(1:21,index1),'b:','linewidth',2), ...
% hold on, plot(HO,x1(1:21,index2),'b:','linewidth',2)
lowval = floor(min(x1(:,index1)));
highval = ceil(max(x1(:,index2)));
axis([0 20 lowval highval])
set(gca,'XTick',0:5:20)
set(gca,'YTick',lowval:0.5:highval)
set(gca,'Layer','top')
title(title_a,'fontsize',10)