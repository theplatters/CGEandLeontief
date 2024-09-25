%% Run Carvalho Gabaix to get sectoral TFP measures 
clear
%cd('/Users/joko/Desktop/Replication Files/GDP Simulatin -- 88 Sector')
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');

%% cREAD JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
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

%% Simulation 
cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:43:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));

%%

A = ones(N,1);
epsilon = .3; % Elasticity of substitution between VA and intermediates
theta = 0.0001; % Elasticity of substitution between intermediates
sigma = 0.4; % Elasticity of substitution between in consumption

epsilon = .6; %
theta = 0.2; % 
sigma = .9; % 

%%
trials = 10;
clear Shocks;
GDP = zeros(trials,1);
parfor_progress(trials);
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha;
init = [ones(N,1)];
%generate idiosyncratic shocks from empirical CDF
%  for j = 1:76
%                 %Shocks(j,:) = exp((emprand((stfp(j,:)-mu(j)),trials,1)));
%                 Shocks(j,:) = exp((randsample((stfp(j,:)-mu(j)),trials,1)));
%  end
%Shocks = datasample(stfp'-repmat(mu,size(stfp',1),1),trials)';
tic
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));

parfor k = 1:trials
            %A = exp(Shocks(:,k));
            %A = Shocks(:,k);
            %A = ones(N,1);
            %A = exp(stfp(:,k)-mu');
            %A = exp(mvnrnd(-1/2*diag(Sigma),Sigma))';         
            A = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))';  
            A = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))';  
            display(A)
            %A = exp(mvnrnd(-diag(Sigma),2*diag(diag(Sigma))))'; % twice as big shocks
            %A(1) = grid(k);
            [Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs_realloc(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),optionsf);         
            display(Soln)

            if exitfl == 1 || exitfl == 2
                GDP(k) = (beta'*Soln.^(1-sigma))^(1/(sigma-1));
                %lambda_simul(:,k) = Soln(1:N)*Soln(1+N,:)/GDP(k);
            end
            parfor_progress;
     
end 
% toc
parfor_progress(0);

correct = (GDP~=0 & log(GDP)>-.4 & log(GDP) <.3); % Only use GDP values that are not caused by numerical errors
mean(log(GDP(correct)));
std(log(GDP(correct)));
skewness(log(GDP(correct)));
kurtosis(log(GDP(correct)))-3;


%save('GDP_histogram_norealloc_normalCOV.mat', 'GDP')
%hist(GDP(GDP~=0),25)
%hist(GDP(abs(GDP)>.1& abs(GDP)<1.4),25)
%% Hulten
Tin = year - 1960;
mean((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))
std((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))
skewness((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))