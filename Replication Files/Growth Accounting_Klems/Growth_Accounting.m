clear
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\Growth Accounting_Klems')
[~, ~, raw] = xlsread('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\Growth Accounting_Klems\usa_wk_mar_2017.xlsx','KLEMdata','A3:J4422');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
year = reshape(data(:,1),68,65);
industry = reshape(data(:,2),68,65);
grossoutput = reshape(data(:,3),68,65);
capital = reshape(data(:,4),68,65);
labor = reshape(data(:,5),68,65);
intermediate = reshape(data(:,6),68,65);
grossoutput1 = reshape(data(:,7),68,65);
capital1 = reshape(data(:,8),68,65);
labor1 = reshape(data(:,9),68,65);
intermediate1 = reshape(data(:,10),68,65);


% Clear temporary variables
clearvars data raw;

% Nonlinear Growth Accounting

labor_share =  (labor./grossoutput);
capital_share =(capital./grossoutput);
intermediate = grossoutput - labor - capital;
inter_share = (intermediate./grossoutput);
GDP = sum(capital+labor,2);
Domar = bsxfun(@rdivide, grossoutput, GDP);
sTFP = (diff(log(grossoutput1))') - labor_share(1:end-1,:)'.*(diff(log(labor1))')...
    - capital_share(1:end-1,:)'.*(diff(log(capital1))')...
    -inter_share(1:end-1,:)'.*(diff(log(intermediate1))');

sVATFP = (diff(log(grossoutput1))') - (labor_share(1:end-1,:)./(labor_share(1:end-1,:)+capital_share(1:end-1,:)))'.*(diff(log(labor1))')...
    - (capital_share(1:end-1,:)./(labor_share(1:end-1,:)+capital_share(1:end-1,:)))'.*(diff(log(capital1))');
    %-inter_share(1:end-1,:)'.*(diff(log(intermediate1))');
VAshares = bsxfun(@rdivide, labor(1:end-1,:)+capital(1:end-1,:),(GDP(1:end-1)));
VATFP = sum(VAshares'.*sVATFP);
TFP = sum(Domar(1:end-1,:)'.*sTFP,1);
Tornqvist = .5*(Domar(1,:)+Domar(end,:));
T = 1;

First_Order_TFP = sum((repmat((Domar(T,:)'),1,67)).*sTFP,1);
Second_Order_TFP = sum((repmat((Tornqvist'),1,67)).*sTFP,1);
plot(year(1:end-1,1), [exp(cumsum(TFP));exp(cumsum(First_Order_TFP)); exp(cumsum(Second_Order_TFP))],'LineWidth',2)
legend('Nonlinear', 'First Order','Second Order', 'Location','SouthEast')
axoptions={'scaled ticks = false',...
           'y tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2, set thousands separator={}}',...
           'x tick label style={/pgf/number format/.cd, set thousands separator={}}',...
           'legend style={font=\scriptsize}'};

matlab2tikz('Cost_Disease.tex','extraAxisOptions',axoptions)
%%
count =1; 
error = zeros(length(TFP)-1,1); 
for step_size = 1:length(TFP)-1
clearvars Hulten_Approx
for T= step_size+1: length(TFP)
      Hulten_Approx(T-step_size)  = Domar(T-step_size,:)*sum(sTFP(:,T-step_size:T-1),2);
end
aggTFP = exp(cumsum(TFP));
%Actual = diff(log(aggTFP([1:step_size:67])));
Actual = [cumsum(TFP')-lagmatrix(cumsum(TFP'),step_size)];
Actual(isnan(Actual)) = [];
error(step_size) = mean(Hulten_Approx)-mean(Actual);
end
%%
for T = 1:step_size:67
    Hulten_4year(count,:) = sum((repmat((Domar(T,:)'),1,67)).*sTFP,1);
    temp(count) = sum(Hulten_4year(count,T:min(T+step_size,67)));
    count = count+1;
end
aggTFP = exp(cumsum(TFP));
mean(-(temp(1:length(Actual))-Actual))
%%
figure; 
for T = 1:67
record(T,:) =  [mean(TFP(1:T)) mean(First_Order_TFP(1:T))];
end
plot([1948:2014]', record);legend('Nonlinear','First Order')
title('Average Growth Rate')
figure
plot([1948:2014]', [TFP' First_Order_TFP']);legend('Nonlinear','First Order')