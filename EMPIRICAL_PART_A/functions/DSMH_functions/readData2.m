nlags = 6;  % this is number of lags in VAR
kcum = 0;   % 0: not cumulated; 1: cumulated for IRFs and var decomps
load data20.txt;   % 1970:Q1 to 2008:Q4
data2=data20;
% col 1 = real GDP 
% col 2 = potential real GDP 
% col 3 = GDP deflator 
% col 4 = CPI 
% col 5 = PCE deflator
% col 6 = fed funds rate (average of 3 months during quarter)
% col 7 = Moody's seasoned Baa corporate bond yield - 10-year Treasury
%         yield
% col 8 = CRB commodity spot price index
% col 9 = average hourly earnings
% col 10 = Dow Jones spot price
% col 11 = real PCE consumption
% col 12 = real investment
% col 13 = IP
% col 14 = unemployment rate
% col 15 = real M2 money stock
% col 16 = trade-weighted US dollar exchange rate: index (broad)
% col 17 = commercial loans
% col 18 = U Mich sentiment index
% col 19 = U Mich inflation expectations
% col 20 = WTI spot price
% col 21 = housing starts
gap = 100*(log(data2(5:end,1)) - log(data2(5:end,2)));
infl1 = 100*(log(data2(5:end,3)) - log(data2(1:end-4,3)));
infl2 = 100*(log(data2(5:end,4)) - log(data2(1:end-4,4)));
infl3 = 100*(log(data2(5:end,5)) - log(data2(1:end-4,5))); % inflation is year over year change in PCE deflator
ffr = data2(5:end,6);
spread = data2(5:end,7);
comm1 = 100*(log(data2(5:end,8)) - log(data2(1:end-4,8)));
wage = 100*(log(data2(5:end,9)) - log(data2(1:end-4,9)));
comm2 = 100*(log(data2(5:end,10)) - log(data2(1:end-4,10)));
cons = 100*(log(data2(5:end,11)) - log(data2(1:end-4,11)));
invest = 100*(log(data2(5:end,12)) - log(data2(1:end-4,12)));
ip = 100*(log(data2(5:end,13)) - log(data2(1:end-4,13)));
unemp = data2(5:end,14);
m2 = 100*(log(data2(5:end,15)) - log(data2(1:end-4,15)));
exch = 100*(log(data2(5:end,16)) - log(data2(1:end-4,16)));
loan = 100*(log(data2(5:end,17)) - log(data2(1:end-4,17)));
senti = data2(5:end,18);
sentid = senti - mean(senti);
inflexp = data2(5:end,19);
oil = 100*(log(data2(5:end,20)) - log(data2(1:end-4,20)));
house = 100*(log(data2(5:end,21)) - log(data2(1:end-4,21)));


   varnames = {' y'; ' \pi'; ' r'};
   shocknames = {' demand'; ' supply'; ' monetary policy'};
   yall = [gap infl3 ffr spread comm1 wage infl1 infl2 comm2 ...
           cons invest ip unemp m2];


tstart = 61;                % start estimation with 1986:Q1
tend = size(yall,1) - 1;    % end estimation with 2008:Q3
YY = yall(tstart:tend,:);
date_start = 1986;
date_end = 2008.5;
time = (1986:0.25:2008.5)';

