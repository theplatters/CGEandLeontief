%% Performs simulation on second order log approximation function of benchmark model. 
%% Run Carvalho Gabaix to get sectoral TFP measures 
clear
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');

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
year = 1982;
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

%% Approximate Second order approximation for benchmark (no-reallocation) calibration
A = ones(N,1);
epsilon = .5; % Elasticity of substitution between VA and intermediates
theta = 0.001; % Elasticity of substitution between intermediates
sigma = 0.9; % Elasticity of substitution between in consumption
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha;
init = [exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A));(beta'*inv(eye(N)-diag(1-alpha)*Omega))'./exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A))];  % judicious choice of starting values
Id = eye(N);
h = 0.000001;
Hess = zeros(N,N);
parfor i = 1:N
    x = A + Id(:,i)*h;  
    Hess(:,i) = (GDP_norealloc_function(x,init, beta, Omega, alpha, epsilon, theta, sigma,L)-GDP_norealloc_function(A,init, beta, Omega, alpha, epsilon, theta, sigma,L))/h;
end
lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
logGDP_2nd = @(X) log(1) + lambda'*log(X) + 1/2*log(X)'*Hess*log(X);

%% Simulation
trials = 50000;
clear Shocks;
GDP = zeros(trials,1);
parfor_progress(trials);
init = [ones(N,1)];

cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:4:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));

parfor k = 1:trials
            %A = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))';  
            A = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))';  
            logGDP(k) = logGDP_2nd(A);
            parfor_progress;
     
end 
toc
parfor_progress(0);



%% Approximate Second order approximation for (full-reallocation) calibration [DIFFERENT METHOD]
A = ones(N,1);
epsilon = .5; % Elasticity of substitution between VA and intermediates
theta = 0.001; % Elasticity of substitution between intermediates
sigma = 0.9; % Elasticity of substitution between in consumption
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha;
init = exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A));  % judicious choice of starting values
Id = eye(N);
h = 0.00001;
parfor_progress(N);
Hess = zeros(N,N);
lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
parfor i = 1:N
    x = A + Id(:,i)*h;  
    p = GDP_realloc_function(x,init, beta, Omega, alpha, epsilon, theta, sigma,L);
    C = (beta'*(p.^(1-sigma)))^(1/(sigma-1));
    q = (Omega*(p.^(1-theta))).^(1/(1-theta));
    y = ((beta'*diag(p)^(-sigma)*C)*inv(eye(N)-diag(p)^epsilon*diag(x)^(epsilon-1)*diag(q)^(theta-epsilon)*diag(1-alpha)*Omega*diag(p)^(-theta)))'; 
    %lambda(i) = (GDP_realloc_function(x,init, beta, Omega, alpha, epsilon, theta, sigma,L)-GDP_realloc_function(TFP,init, beta, Omega, alpha, epsilon, theta, sigma,L))/delta;
    lambda_new = p.*y/C;
    Hess(:,i) = (lambda_new-lambda)/h;
     parfor_progress;
end

logGDP_2nd = @(X) log(1) + lambda'*log(X) + 1/2*log(X)'*Hess*log(X);
%% Simulation for full reallocation case
trials = 50000;
clear Shocks;
GDP = zeros(trials,1);
cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:4:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));
parfor_progress(trials);
init = [ones(N,1)];
parfor k = 1:trials
            %A = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))';  
            A = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))';  
            logGDP(k) = logGDP_2nd(A);
            parfor_progress;
     
end 
toc
parfor_progress(0);

mean(logGDP)
std(logGDP)
skewness(logGDP)
kurtosis(logGDP)-3

