%% Runs simulations with adjustment cost in production and consumption. 
% step_size determines size of shocks. Full reallocation of labor. 
%% Run Carvalho Gabaix to get sectoral TFP measures 
 
clear
 
%cd('C:\Users\Maria\Documents\Baqaee files')
run('C:\Users\BAQAEE\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp aggtfp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');


for round = 1:4
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')

step_size_grid = [1;4;1;4];
theta_grid = [0.001;0.001;0.2;0.2];
epsilon_grid = [0.5; 0.5; 0.6; 0.6];
kappa_grid = [3;3;4;4];

step_size = step_size_grid(round); 

cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:step_size:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));

 
%% READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix

cd('C:\Users\BAQAEE\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')
load us80dbasedata.mat
cd('C:\Users\BAQAEE\Dropbox\Work\Hultens Theorem\Calibration\Stuck Intermediates and Adjustment Costs')
 
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
alpha = (vadd(:,year-1959)./grossy(:,year-1959)); % set the factor share by industry
 
N = length(Omega);
N_orig = N; 
 
beta = grossy(:,year-1959)'*(eye(N)-diag(1-alpha)*Omega);
beta(beta<0) = 0; %remove industries with negative implied final sales
beta = beta/sum(beta); %normalize consumption vector to sum to unity. 
beta = beta';
 
 
%% Simulation with no reallocation
A = ones(N,1);
epsilon = epsilon_grid(round)*ones(N, 1); % Elasticity of substitution between VA and intermediates (NOTE: epsilon and theta reversed relative to notation in the paper)
theta = theta_grid(round)*ones(N, 1); % Elasticity of substitution between intermediatestoc
sigma = .9; % Elasticity of substitution between in consumption
eta  = ones(N, 1);
 
%% Fmincon Options
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
 
%% Relabel # 1 -- computing the endoments
M = 1;

[Omega_tr, Elasticities_tr, beta_tr] = Relabel_HT_adj_costs_realloc(Omega, alpha, epsilon, theta, beta, N, M, sigma);
N = length(beta_tr);
b = beta_tr;
%Omega_tr = bsxfun(@rdivide,Omega_tr, sum(Omega_tr,2));
Omega_tr(isnan(Omega_tr))=0;
epsilon_tr = Elasticities_tr;
 
alpha_tr = zeros(N,1);
alpha_tr(end-M+1:end, :) = 1;
A_tr = ones(N, 1);
 
L_inv  = inv(eye(N)-Omega_tr);
Matr = (beta_tr'*L_inv)';
 
L = Matr.*alpha_tr; % steady-state allocation of labor 
L = L(end-M+1:end);
shares = Matr.*ones(N, 1);
 
init2 = shares([1:3*N_orig+2, end-1:end]); 
 
%% Simulation

%rng(256)
trials = 10000; % number of draws
clear Shocks LShocks;
GDP = zeros(trials,1);
lost_res = zeros(trials, 1);
profit_share = zeros(trials,1);

L_inv  = inv(eye(N)-Omega_tr);
L_inv = sparse(L_inv);
 
L_orig = init2(end-M+1:end); % steady-state allocation of labor 
init_orig = [ones(N,1);(beta'*inv(eye(N_orig)-diag(1-alpha)*Omega))']; 
 
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
 
bar_X = (1 - alpha).*shares(1+2:N_orig+2);
kappa = kappa_grid(round)*ones(N_orig+1, 1);
% Specify the adjustment costs for consumption here (if different from production)
%kappa(1) = 5; %adjustment costs on consumption goods. 

parfor_progress(trials);
parfor k = 1:trials
            k
            %A_orig = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))'; % independent log normal shocks
            A_orig = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))'; 
            A = vertcat(1, 1, A_orig, ones(2*N_orig+M,1));
            
            init_1 = (exp(-L_inv*log(A)));
            init_2 =  (init2);
            init = vertcat(init_1(1:N-M), init_2(1:N-M));  % judicious choice of starting values
            
            tic
            %[Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Derivs_HT_new_reallocI(X, A, beta_tr, Omega_tr, epsilon_tr, 1, L, N, N_orig, M, bar_X, kappa), optionsf);
            [Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Derivs_HT_new_reallocI(X, A, beta_tr, Omega_tr, epsilon_tr, 1, L, N, N_orig, M, bar_X, kappa),zeros(length(init),1),[],'Knitro_options.opt');
            toc

            %if exitfl == 1 || exitfl == 2% if the solver didn't give error
            if exitfl == 0% solver no error (knitro)                                                  
                w  = 1;
                
                p = Soln(1:N-M);
                y = Soln(N-M+1:2*(N-M));
                
                GDP(k) = y(1);
                lost_res(k) = (p(1)*(y(2)-y(1))+p(N_orig+1+2:2+2*N_orig)'*(y(2+2*N_orig+1:2+3*N_orig)-y(2+N_orig+1:2+2*N_orig)))/(p(1)*y(1));
                profit_share(k) = (y(1)*p(1)-y(2)*p(2)+p(N_orig+1+2:2+2*N_orig)'*y(2+N_orig+1:2+2*N_orig) - p(2+2*N_orig+1:2+3*N_orig)'*y(2+2*N_orig+1:2+3*N_orig))/(p(1)*y(1));
                lambda_simul(:,k) = p.*y/(p(1)*y(1));
            end
            parfor_progress;
end 
parfor_progress(0); 
 
correct = (GDP~=0 & log(GDP)>-.4 & log(GDP) <.3); % Only use GDP values that are not caused by numerical errors

mean_GDP = mean(log(GDP(correct)));
std_GDP = std(log(GDP(correct)));
skw_GDP = skewness(log(GDP(correct)));
kurt_GDP = kurtosis(log(GDP(correct)))-3;

domar_weights = shares(1+2:N_orig+2);
 

lambd = lambda_simul(1+2:N_orig+2, correct);

volatility = sqrt(var(lambd, 0, 2)')*domar_weights(1:N_orig)./sum(domar_weights);

domar_vol = domar_weights'*std(bsxfun(@minus, log(lambd)',log(domar_weights)'))'

cd('Results with Kappa')
save(sprintf('%s%s%g%s%g%s%g%s%g%s%g%s', 'Realloc_GDP','sigma',sigma,'epsilon', epsilon(1),'theta', theta(1), 'step_size',step_size,'kappa',kappa(2),'.mat'))
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\Stuck Intermediates and Adjustment Costs')
end

%% Runs simulations with adjustment costs for consumption and production. 
% Labor is non-reallocatable.

%% Run Carvalho Gabaix to get sectoral TFP measures 
 
clear
 


%cd('C:\Users\Maria\Documents\Baqaee files')
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp aggtfp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');

for step_size = 1:3:4
cd('C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector')

    
cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:step_size:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));


%% READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix
 
load us80dbasedata.mat
 cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\Stuck Intermediates and Adjustment Costs')
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
alpha = (vadd(:,year-1959)./grossy(:,year-1959)); % set the factor share by industry
 
N = length(Omega);
N_orig = N; 
 
beta = grossy(:,year-1959)'*(eye(N)-diag(1-alpha)*Omega);
beta(beta<0) = 0; %remove industries with negative implied final sales
beta = beta/sum(beta); %normalize consumption vector to sum to unity. 
beta = beta';
 
 
%% Simulation with no reallocation
A = ones(N,1);
epsilon = .6*ones(N, 1); % Elasticity of substitution between VA and intermediates (NOTE: epsilon and theta reversed relative to notation in the paper)
theta = 0.2*ones(N, 1); % Elasticity of substitution between intermediatestoc
sigma = .9; % Elasticity of substitution between in consumption
eta  = ones(N, 1);
 
%% Fmincon Options
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
 
%% Relabel # 1 -- computing the endoments
M = N_orig;

[Omega_tr, Elasticities_tr, beta_tr] = Relabel_HT_adj_costsI(Omega, alpha, epsilon, theta, beta, N, M, sigma);
N = length(beta_tr);
b = beta_tr;
%Omega_tr = bsxfun(@rdivide,Omega_tr, sum(Omega_tr,2));
Omega_tr(isnan(Omega_tr))=0;
epsilon_tr = Elasticities_tr;
 
alpha_tr = zeros(N,1);
%reset alphas for full/no reallocation
alpha_tr(end-M+1:end, :) = 1;
A_tr = ones(N, 1);
 
L_inv  = inv(eye(N)-Omega_tr);

Matr = (beta_tr'*L_inv)';
 
L = Matr.*alpha_tr; % steady-state allocation of labor 
L = L(end-M+1:end);
shares = Matr.*ones(N, 1);
 
init2 = shares([1:3*N_orig+2, end-N_orig+1:end]); 
 
%% Simulation

%rng(256)
trials = 10000; % number of draws
clear Shocks LShocks;
GDP = zeros(trials,1);
lost_res = zeros(trials, 1);
profit_share = zeros(trials, 1);

L_inv  = inv(eye(N)-Omega_tr);
L_inv = sparse(L_inv);
 
L_orig = init2(end-M+1:end); % steady-state allocation of labor 
init_orig = [ones(N,1);(beta'*inv(eye(N_orig)-diag(1-alpha)*Omega))']; 
 
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
 parfor_progress(trials);
bar_X = (1 - alpha).*shares(1+2:M+2);

% Adjustment costs for all sectors are the same
kappa = 2*ones(M, 1);

parfor k = 1:trials
            k
            A_orig = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))'; 
            A = vertcat(1, 1, A_orig, ones(3*N_orig,1));
            
            init_1 = (exp(-L_inv*log(A)));
            init_2 =  (init2);
            init = vertcat(init_1(1:N-M), init_2(1:N-M));  % judicious choice of starting values
            %tic
            %[Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Derivs_HT_newVII(X, A, beta_tr, Omega_tr, epsilon_tr, 1, L, N, M, bar_X, kappa), optionsf);
            [Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Derivs_HT_newVII(X, A, beta_tr, Omega_tr, epsilon_tr, 1, L, N, M, bar_X,kappa),[],[], 'Knitro_options.opt');
            %toc

            %if exitfl == 1 || exitfl == 2% if the solver didn't give error
            if exitfl == 0% solver no error (knitro)                                 
                 
                w  =  ((Omega_tr(1:N-M, N-M+1:N)'.*(1./repmat(L, 1, N-M))).^(repmat(1./epsilon_tr(1:N-M)', M, 1)))*(Soln(1+N-M:2*(N-M)).^(1./epsilon_tr(1:N-M)).*Soln(1:(N-M)).*(A(1:N-M).^((epsilon_tr(1:N-M)-1)./epsilon_tr(1:N-M))));
                w(isnan(w)) = 0;
                
                p = Soln(1:N-M);
                y = Soln(N-M+1:2*(N-M));
                
                GDP(k) = y(1);
                lost_res(k) = (p(1)*(y(2)-y(1))+p(N_orig+1+2:2+2*N_orig)'*(y(2+2*N_orig+1:2+3*N_orig)-y(2+N_orig+1:2+2*N_orig)))/(p(1)*y(1));
                profit_share(k) = (y(1)*p(1)-y(2)*p(2)+p(N_orig+1+2:2+2*N_orig)'*y(2+N_orig+1:2+2*N_orig) - p(2+2*N_orig+1:2+3*N_orig)'*y(2+2*N_orig+1:2+3*N_orig))/(p(1)*y(1));
                lambda_simul(:,k) = p.*y/(p(1)*y(1));
            end
          parfor_progress;
end 
 parfor_progress(0);
 
correct = (GDP~=0 & log(GDP)>-.4 & log(GDP) <.4); % Only use GDP values that are not caused by numerical errors

mean_GDP = mean(log(GDP(correct)));
std_GDP = std(log(GDP(correct)));
skw_GDP = skewness(log(GDP(correct)));
kurt_GDP = kurtosis(log(GDP(correct)))-3;

domar_weights = shares(1+2:N_orig+2);
 
%for i = 1:trials
%    domar_w(1:N_orig, i) = sum(diag(lambda_simul(N_orig+1:2*N_orig, i))*(Omega*diag(1- zeta)), 1)';
%end

lambd = lambda_simul(1+2:N_orig+2, correct);

%volatility = sqrt(var(lambd, 0, 2)')*domar_weights(1:N_orig)./sum(domar_weights);
    
domar_vol = domar_weights'*std(bsxfun(@minus, log(lambd)',log(domar_weights)'))'

cd('Results with Kappa')
save(sprintf('%s%s%g%s%g%s%g%s%g%s%g%s', 'NoRealloc_GDP','sigma',sigma,'epsilon', epsilon(1),'theta', theta(1), 'step_size',step_size,'kappa',kappa(2),'.mat'))
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\Stuck Intermediates and Adjustment Costs')
end