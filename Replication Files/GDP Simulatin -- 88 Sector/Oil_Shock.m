 %% Run Carvalho Gabaix to get sectoral TFP measures 
clear
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp aggtfp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');
stfp_fe =  csvread('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\stfp_resids_fe.csv', 1,0);
stfp_tt =  csvread('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\stfp_resids_linear_trend.csv', 1,0);
stfp_tt_ind =  csvread('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\stfp_resids_ttind.csv', 1,0);
%% READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix


load us80dbasedata.mat

%AUXILIARY COUNTERS TO TURN ORIGINAL DATA INTO MATRIX(SECTORxYEAR)
startcount=1:46:4048;
for k=1:size(startcount,2)
    endcount(k)=startcount(k)-1;
end
endcount(89)=size(data,1);

%NOMINAL GROSS OUTPUT IS THE SUM OF NOMINAL CAPITAL, LABOR AND ALL
%INTERMEDIATE INPUTS (INCLUDING NONCOMPETING IMPORTS)
for i=1:size(startcount,2)
    grossy(i,:)=data(startcount(i):endcount(i+1),3)';
end

%NOMINAL VALUE ADDED IS THE SUM OF NOMINAL CAPITAL AND LABOR
for i=1:size(startcount,2)
    vadd(i,:)=(data(startcount(i):endcount(i+1),4)')+(data(startcount(i):endcount(i+1),5)');
end

%NOMINAL CAPITAL
for i=1:size(startcount,2)
    capital(i,:)=(data(startcount(i):endcount(i+1),4)');
end

%NOMINAL LABOR
for i=1:size(startcount,2)
    labor(i,:)=(data(startcount(i):endcount(i+1),5)');
end

%REMOVE GOVERNMENT SECTORS & RENTS IMPUTED FROM OWNER-OCCUPIED HOUSING
temp=80:1:88;
temp=[60;temp'];
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];
clear temp

%REMOVE SECTORS FOR WHICH THERE IS NO GROSS SALES (8, URANIUM ORES AND 62,
%RENTING OF MACHINERY)

[temp]=find(sum(grossy')==0);
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];

%% Initialize a year for IO Matrix
year = 1971;
temp=80:1:88;
temp=[8;60;62;temp'];
IO = data(find(data(:,1)==year),:);
IO(:,[1 3 4 5 94]) = []; % delete year, gross output, capital, labor, noncompetitive imports
IO(temp,:) = []; % remove government sectors, and sectors with no gross sales
Ind = IO(:,1); %store industry names
IO(:,1) = [];
IO(:,temp) =[];
Omega = bsxfun(@rdivide, IO, sum(IO,2));
%Omega = diag(1-vadd(:,year-1959)./grossy(:,year-1959))*Omega; % scale IO table by intermediate input share
alpha = (vadd(:,year-1959)./grossy(:,year-1959)); % set the factor share by industry
N = length(Omega);
beta = grossy(:,year-1959)'*(eye(N)-diag(1-alpha)*Omega);
beta(beta<0) = 0; %remove industries with negative implied final sales
beta = beta/sum(beta); %normalize consumption vector to sum to unity. 
beta = beta';


%% Shocks to Oil Industry
Ind = 7; % Industry with most action is 7, least action is 17, sames sales as oil is 53;
A = ones(N,1);
epsilon = .02; % Elasticity of substitution between VA and intermediates
theta = 0.25; % Elasticity of substitution between intermediates
sigma = .25; % Elasticity of substitution between in consumption
lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
clear Shocks LShocks GDP;
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha;
init = [ones(N,1);(beta'*inv(eye(N)-diag(1-alpha)*Omega))'];
%generate idiosyncratic shocks from empirical CDF
tic
A = ones(N,1);    
A(Ind) = 0.8*1/2+1/2;
init = [exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A));(beta'*inv(eye(N)-diag(1-alpha)*Omega))'./exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A))];  % judicious choice of starting values
[Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),[],[],'Knitro_options.opt');
if exitfl == 0
          GDP_oil = (Soln(1:N).*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(Soln(N+1:2*N).^(1/epsilon)).*(1./L).^(1/epsilon))'*L;
            lambda_simul = Soln(1:N)*Soln(1+N,:)/GDP_oil;
end
