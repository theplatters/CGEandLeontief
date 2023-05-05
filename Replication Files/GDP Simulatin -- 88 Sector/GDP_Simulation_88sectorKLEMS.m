%% Run Carvalho Gabaix to get sectoral TFP measures
clear
%cd('/Users/Joko/Desktop/Replication Files/GDP Simulatin -- 88 Sector')
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp aggtfp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');
%% READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix


load us80dbasedata.mat


grossy = reshape(data(:,3),46,88)'; %NOMINAL OUTPUT
capital = reshape(data(:,4),46,88)'; %NOMINAL CAPITAL
labor = reshape(data(:,5),46,88)'; %NOMINAL LABOR
vadd = labor + capital; %NOMINAL VALUE



%REMOVE GOVERNMENT SECTORS & RENTS IMPUTED FROM OWNER-OCCUPIED HOUSING

temp=[60 ,80:88];
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];
clear temp

%%

%REMOVE SECTORS FOR WHICH THERE IS NO GROSS SALES (8, URANIUM ORES AND 62,
%RENTING OF MACHINERY)

[temp]=find(sum(grossy')==0);
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];

%% Initialize a year for IO Matrix
year = 1982;
temp=80:88;
temp=[8;60;62;temp'];
IO = data(find(data(:,1)==year),:);
IO(:,[1 3 4 5 94]) = []; % delete year, gross output, capital, labor, noncompetitive imports
IO(temp,:) = []; % remove government sectors, and sectors with no gross sales
Ind = IO(:,1); %store industry names
IO(:,1) = [];
IO(:,temp) =[];

Omega = bsxfun(@rdivide, IO, sum(IO,2)); %Divides each elememt IO[i,j] with sum_i IO[i,j]
alpha = (vadd(:,year-1959)./grossy(:,year-1959)); %set the factor share by industry
N = length(Omega);


beta = grossy(:,year-1959)'*(eye(N)-diag(1-alpha)*Omega);
beta(beta<0) = 0; %remove industries with negative implied final sales
beta = beta/sum(beta); %normalize consumption vector to sum to unity.
beta = beta';


domar_weights = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';

%% Fmincon Options
optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
%optionsf = optimoptions('Algorithm','interior-point')

%% Simulation with no reallocation
% create 4-yearly standard deviation of shocks
cum_stfp = cumsum(log(1+stfp)')';
cum_stfp_4year = cum_stfp(:,[1:4:size(stfp,2)]);
Sigma_4year = cov(diff(cum_stfp_4year'));
%%




A = ones(N,1);

% WHAT
%epsilon = .3; % Elasticity of substitution between VA and intermediates (NOTE: epsilon and theta reversed relative to notation in the paper)
%epsilon = .15; % To match volatilities
%theta = 0.0001; % Elasticity of substitution between intermediates
%sigma = .4; % Elasticity of substitution between in consumption

epsilon = .5; %
theta = 0.001; %
sigma = .9; %


trials = 1; % number of draws
clear Shocks LShocks;
GDP = zeros(trials,1);
variances = (movingvar(stfp',5)'); % rolling estimate of variance of TFP
var_cri = variances(:,22);

var_cri = diag(Sigma)*2; % crisis episode variances
%parfor_progress(trials);
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha; % steady-state allocation of labor
init = [ones(N,1);(beta'*inv(eye(N)-diag(1-alpha)*Omega))']; % initial conditions for the solver

tic

% Oil Shocks -- these aren't used anymore
Oil = ones(45,1);
Oil(13:15) = exp(stfp(7,13:15)-mu(7));
Oil(19:21) = exp(stfp(7,19:21)-mu(7));


parfor k = 1:trials
            A = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma))))';
            % A = exp(mvnrnd(-1/2*diag(Sigma_4year),diag(diag(Sigma_4year))))';
             init = [exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A));(beta'*inv(eye(N)-diag(1-alpha)*Omega))'./exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(A))];  % judicious choice of starting values
            [Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),optionsf);


            if exitfl == 1 || exitfl == 2% solver no error (fmincon)
                GDP(k) = (Soln(1:N).*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(Soln(N+1:2*N).^(1/epsilon)).*(1./L).^(1/epsilon))'*L;
                lambda_simul(:,k) = Soln(1:N).*Soln(1+N:end)/GDP(k);
            end
            %parfor_progress;

end
toc


parfor_progress(0);


correct = (GDP~=0 & log(GDP)>-.4 & log(GDP) <.3); % Only use GDP values that are not caused by numerical errors
mean(log(GDP(correct)));
std(log(GDP(correct)));
skewness(log(GDP(correct)));
kurtosis(log(GDP(correct)))-3;



lambd = lambda_simul(:, correct);

volatility = sqrt(var(lambd, 0, 2)')*domar_weights(:)./sum(domar_weights);



% record domar weightrs over time
count = 1;
for year = 1960:2005
%year = 1960;
  temp=[8,60,62,80:88];

  IO = data(find(data(:,1)==year),:);
  IO(:,[1 3 4 5 94]) = []; % delete year, gross output, capital, labor, noncompetitive imports
  IO(temp,:) = []; % reove government sectors, and sectors with no gross sales
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
  lambda(:,count) = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
  count = count+1;
end

mean(lambda')*(std(diff(log(lambda_simul(:,correct))')))'


%% Hulten
Tin = year - 1960;
Shocks = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma)),trials))';
mean((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))
std((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))
skewness((grossy(:,Tin)./sum(vadd(:,Tin)))'*log(Shocks))

lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
trials = 50000;
Shocks = exp(mvnrnd(-1/2*diag(Sigma),diag(diag(Sigma)),trials))';
GDP_CD = exp(lambda'*log(Shocks))';
save('GDP_simulation_50K_CD.mat', 'GDP_CD')

%% Shocks to single industries and comparison to Hulten
Ind = 7; % Industry with most action is 7, least action is 17, sames sales as oil is 53;
A = ones(N,1);
a = 0.7;
b = 1.3;
lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
epsilon = .3; % Elasticity of substitution between VA and intermediates
theta = 0.0001; % Elasticity of substitution between intermediates
sigma = .4; % Elasticity of substitution between in consumption

%
 epsilon = .5; % Elasticity of substitution between VA and intermediates
 theta = 0.0001; % Elasticity of substitution between intermediates
 sigma = .9; % Elasticity of substitution between in consumption


clear Shocks LShocks GDP;
L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha;
init = [ones(N,1);(beta'*inv(eye(N)-diag(1-alpha)*Omega))'];
%generate idiosyncratic shocks from empirical CDF
tic
M = 10;
GDP = zeros(2*M,1);
grid = linspace(1.0, a, M);
list = [7;53;8]
for k = 1:length(list)
      Ind = list(k);
      grid = linspace(1.0, a, M);
  for j = 1:M
      A = ones(N,1);
      A(Ind) = grid(j);
      [Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),[],[],'Knitro_options.opt');
      if exitfl == 0
        GDP(j,k) = (Soln(1:N).*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(Soln(N+1:2*N).^(1/epsilon)).*(1./L).^(1/epsilon))'*L;
      end
      init = Soln;
  end


  grid = linspace(1.0, b, M);
  init = [ones(N,1);(beta'*inv(eye(N)-diag(1-alpha)*Omega))'];
  for j = 1:M
          A = ones(N,1);
          A(Ind) = grid(j);
              [Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs(X,  A, beta, Omega, alpha, epsilon, theta, sigma,L),[],[],'Knitro_options.opt');
              if exitfl == 0
                  GDP(j+M,k) = (Soln(1:N).*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(Soln(N+1:2*N).^(1/epsilon)).*(1./L).^(1/epsilon))'*L;
              end
              init = Soln;
  end
  GDP(1:M,k) = flipud(GDP(1:M,k));
end
toc
GDP(M,:) =[];

figure;
plot(linspace(a,b,2*M-1), [log(GDP(:,1)) lambda(list(1))*(log((linspace(a,b,2*M-1))))'],'LineWidth', 2)
legend('Model', 'Hulten','Location', 'SouthEast')
xlabel('log(TFP)')
ylabel('log(GDP)')
axis([a b -.2 .2])
%% Oil V Retail
close all
f = figure;
a1 = axes('Parent',f);
histogram(stfp(:),'Normalization','probability','BinLimits', [log(a) log(b)])
a2 = axes('Parent',f);
h = zeros(4,1);

h(1) = plot(log(linspace(a,b,2*M-1)), [log(GDP(:,1))],'LineWidth', 2,'Color', 'r'); hold on;
h(2) = plot(log(linspace(a,b,2*M-1)), [log(GDP(:,2))],'LineWidth', 2,'Color','black');
h(3) = plot(log(linspace(a,b,2*M-1)), bsxfun(@times,lambda(list(1)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','r');
h(4) = plot(log(linspace(a,b,2*M-1)), bsxfun(@times,lambda(list(2)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','black');

plot(log(linspace(a,b,2*M-1)), [log(GDP(:,1))],'LineWidth', 2,'Color','r'); hold on;
plot(log(linspace(a,b,2*M-1)), [log(GDP(:,2))],'LineWidth', 2,'Color','black');

plot(log(linspace(a,b,2*M-1)), bsxfun(@times,lambda(list(1)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','r');
plot(log(linspace(a,b,2*M-1)), bsxfun(@times,lambda(list(2)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','black');

legend(h,'Oil','Retail','Hulten Oil','Hulten Retail','Location','NorthWest')
hold off
xlabel('log(TFP)')
ylabel('log(Y/Y)')
set(a2,'Color','none')
set(a2,'YAxisLocation','left')
set(a1,'XTickLabels',[])
set(a1,'YTickLabels',[])
set(a1,'YTick',[])
set(a1,'XTick',[])
set(a1, 'legend','off')
matlab2tikz('OilvRetail.tex')

%% Oil V. Construction

figure;
h = zeros(4,1);
h(1) = plot(linspace(a,b,2*M-1), [log(GDP(:,1))],'LineWidth', 2,'Color', 'black'); hold on;
h(2) = plot(linspace(a,b,2*M-1), [log(GDP(:,3))],'LineWidth', 2,'Color','green');
h(3) = plot(linspace(a,b,2*M-1), bsxfun(@times,lambda(list(1)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','black');
h(4) = plot(linspace(a,b,2*M-1), bsxfun(@times,lambda(list(3)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','green');

plot(linspace(a,b,2*M-1), [log(GDP(:,1))],'LineWidth', 2,'Color','black'); hold on;
plot(linspace(a,b,2*M-1), [log(GDP(:,3))],'LineWidth', 2,'Color','green');
plot(linspace(a,b,2*M-1), bsxfun(@times,lambda(list(1)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','black');
plot(linspace(a,b,2*M-1), bsxfun(@times,lambda(list(3)),(log((linspace(a,b,2*M-1)))))','--','LineWidth', 2,'Color','green');
legend(h,'Oil','Construction','Hulten Oil','Hulten Construction','Location','SouthEast')
hold off
xlabel('TFP')
ylabel('log(GDP)')
matlab2tikz('OilvConstruction.tex')
