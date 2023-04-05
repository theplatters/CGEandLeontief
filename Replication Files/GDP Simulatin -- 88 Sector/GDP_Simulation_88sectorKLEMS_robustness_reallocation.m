%% This code runs simulation over permutations of elasticities with full reallocation
% recroding the results. It can run with 1 year shocks or 4 year shocks. 
%% Run Carvalho Gabaix to get sectoral TFP measures 


run('industries_data_privatesector_1960_2008.m');
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry


theta_grid = [.001;0.2;0.99];
sigma_grid = [.8;0.9;0.99];
epsilon_grid = [.5;0.4;.6;0.99];

step_size = 4; %number of years for size of shocks.

for count1 = 1:length(epsilon_grid)
    for count2 = 1:length(theta_grid)
        for count3 = 1:length(sigma_grid)
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')
clearvars -except stfp agggdp aggtfp count1 count2 count3 theta_grid epsilon_grid sigma_grid step_size

%% READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix

Sigma = cov(stfp');
mu = mean(stfp');
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
clearvars temp

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
alpha = (vadd(:,year-1959)./grossy(:,year-1959)); % set the factor share by industry
N = length(Omega);
beta = grossy(:,year-1959)'*(eye(N)-diag(1-alpha)*Omega);
beta(beta<0) = 0; %remove industries with negative implied final sales
beta = beta/sum(beta); %normalize consumption vector to sum to unity. 
beta = beta';
%% Fmincon Options
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));

%% Simulation with no reallocation
cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:step_size:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));

A = ones(N,1);
%epsilon = .3; % Elasticity of substitution between VA and intermediates (NOTE: epsilon and theta reversed relative to notation in the paper)
epsilon = epsilon_grid(count1); % To match volatilities
theta = theta_grid(count2); % Elasticity of substitution between intermediates
sigma = sigma_grid(count3); % Elasticity of substitution between in consumption
trials = 50000; % number of draws
clearvars Shocks LShocks;
GDP = zeros(trials,1);
variances = (movingvar(stfp(:,:)',5)'); % rolling estimate of variance of TFP
var_cri = variances(:,22);
var_cri = diag(Sigma)*2; % crisis episode variances
parfor_progress(trials);
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha; % steady-state allocation of labor 
init = [ones(N,1);(beta'*inv(eye(N)-diag(1-alpha)*Omega))']; % initial conditions for the solver

% below, draw shocks from the empirical CDF of TFP
 %for j = 1:76
                %LShocks(j,:) = (emprand((stfp(j,:)-mu(j)),trials,1));
 %               Shocks(j,:) = exp((randsample((stfp(j,:)-mu(j)),trials,1)));
 %end
%Shocks = datasample(stfp'-repmat(mu,size(stfp',1),1),trials)';
tic

%Shocks = exp(bsxfun(@minus, LShocks', mean(LShocks')))';
% wage is given by p.*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon)).*(1./L).^(1/epsilon)
% Oil Shocks
Oil = ones(45,1);
Oil(13:15) = exp(stfp(7,13:15)-mu(7));
Oil(19:21) = exp(stfp(7,19:21)-mu(7));

%Sigma = cov(stfp_tt');
parfor k = 1:trials
            %A = exp(Shocks(:,k));
            %A = Shocks(:,k);
            %A = ones(N,1);
            %A = exp(stfp(:,k)-mu');
            %A = exp(stfp_fe(:,k)-mu');
            %A = exp(stfp_tt(:,k)-mu');
            %A = exp(mvnrnd(-1/2*diag(Sigma),Sigma))';   
            %rng(234)
            %A = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))';  
            A = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))';  
            %A = ones(N,1);
            %A(7) = Oil(k);
            %A = exp(mvnrnd(-1/2*var_cri,diag((var_cri))))';  
            %A = exp(mvnrnd(-diag(Sigma),2*diag(diag(Sigma))))'; % twice as big shocks
            %A(1) = grid(k);
            init = [exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A))];  % judicious choice of starting values
            
            % pick solver type
            %[Soln,~,exitfl] = knitromatlab(@(X) trivial(X),ones(N,1),[],[],[],[],[],[], @(X)Simulation_Derivs_realloc(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),[],[],'Knitro_options.opt');
            [Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs_realloc(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L), optionsf);
           
            if exitfl == 1 || exitfl == 2% solver no error (fmincon)
            %if exitfl == 0% solver no error (knitro)
                GDP(k) = (beta'*Soln.^(1-sigma))^(1/(sigma-1));
                prices(:,k) = Soln; 
                q = (Omega*(Soln.^(1-theta))).^(1/(1-theta));
                lambda_simul(:,k) = (beta.*(Soln.^(1-sigma)))'*inv(eye(N)-diag(1-alpha)*diag((q./Soln).^(1-epsilon))*diag(q.^(theta-epsilon))*diag(A.^(epsilon-1))*Omega*diag((Soln).^(1-theta)));
            end
     
          
            parfor_progress;
     
end 
toc
parfor_progress(0);
correct = (GDP~=0 & log(GDP)>-.4 & log(GDP) <.3); % Only use GDP values that are not caused by numerical errors
mean(log(GDP(correct)));
std(log(GDP(correct)));
skewness(log(GDP(correct)));
kurtosis(log(GDP(correct)))-3;


domar_weights = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';

lambd = lambda_simul(:, correct);

volatility = sqrt(var(lambd, 0, 2)')*domar_weights(:)./sum(domar_weights);
cd(sprintf('%s%g%s', 'C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\',step_size,' year\Full Reallocation'))
save(sprintf('%s%s%g%s%g%s%g%s', 'GDP','sigma',sigma,'epsilon', epsilon,'theta', theta,'.mat'))
[count1 count2 count3]
        end
    end
end

