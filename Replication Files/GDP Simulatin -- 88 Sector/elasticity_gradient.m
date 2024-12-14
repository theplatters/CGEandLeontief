%% Run Carvalho Gabaix to get sectoral TFP measures
clear
%cd('/Users/Joko/Desktop/Replication Files/GDP Simulatin -- 88 Sector')
run('industries_data_privatesector_1960_2008.m');
clearvars -except stfp agggdp aggtfp
stfp(:,1) = []; % remove the first year growth which is blank
stfp(end,:) = []; % remove the private household industry
Sigma = cov(stfp');
mu = mean(stfp');

save("stfp.mat","stfp")
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
year = 1980;
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

L = (beta'*inv(eye(N)-diag(1-alpha)*Omega))'.*alpha; % steady-state allocation of labor

epsilons = linspace(0.015,0.99);
thetas  = linspace(0.015,0.99);
sigmas = linspace(0.015,0.99);

sectors = [10, 20, 23];


gdp_eps = zeros(3,100);
gdp_theta = zeros(3,100);
gdp_sigma = zeros(3,100);


mp_eps = zeros(3,100);
mp_theta = zeros(3,100);
mp_sigma = zeros(3,100);

sigma_orig = 0.015;
theta_orig = 0.015;
eps_orig = 0.015;
%%
i = 1;
for sector = sectors
    shock = ones(N,1);
    shock(sector) = 1.1;
    
    epsilon = eps_orig;
    theta = theta_orig;
    j = 1;
    for sigma = sigmas
        [gdp,mp] = eg(shock,epsilon,theta,sigma,beta,Omega,alpha,L,N);
        gdp_sigma(i,j) = gdp;
        mp_sigma(i,j) = mp;

        j = j + 1;
    end

    epsilon = eps_orig;
    sigma = sigma_orig;
    j = 1;
     for theta = thetas
        [gdp,mp] = eg(shock,epsilon,theta,sigma,beta,Omega,alpha,L,N);
        gdp_theta(i,j) = gdp;
        mp_theta(i,j) = mp;

        j = j +1;
     end
    theta = theta_orig;
    sigma = sigma_orig;
    j = 1;
     for epsilon = epsilons
        [gdp,mp] = eg(shock,epsilon,theta,sigma,beta,Omega,alpha,L,N);
        gdp_eps(i,j) = gdp;
        mp_eps(i,j) = mp;
        j = j + 1;
     end
     i = i+1
end


%%
% Plotting
x = linspace(0.015,0.99)
% Plot the lines
h = figure; % Open a new figure window
plot(x, gdp_theta(1,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, gdp_sigma(1,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, gdp_eps(1,1:100), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot
saveas(h,"eg_supply_1_0015",'png')

%%
h = figure; % Open a new figure window
plot(x, mp_theta(1,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, mp_sigma(1,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, mp_eps(1,:), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot

saveas(h,"eg_supply_mp_1_0015",'png')

%%
% Sector 2
x = linspace(0.015,0.99)
% Plot the lines
h = figure; % Open a new figure window
plot(x, gdp_theta(2,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, gdp_sigma(2,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, gdp_eps(2,:), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot

saveas(h,"eg_supply_2_0015",'png')
%%
h = figure; % Open a new figure window
plot(x, mp_theta(2,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, mp_sigma(2,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, mp_eps(2,:), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot
saveas(h,"eg_supply_2_mp_0015",'png')

%%
% Sector 3
x = linspace(0.015,0.99)
% Plot the lines
figure; % Open a new figure window
plot(x, gdp_theta(3,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, gdp_sigma(3,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, gdp_eps(3,:), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot

saveas(h,"eg_supply_3_0015",'png')

%%              
h = figure; % Open a new figure window
plot(x, mp_theta(3,:)  , 'r', 'LineWidth', 1.5); % First line in red
hold on; % Hold the current plot
plot(x, mp_sigma(3,:), 'g', 'LineWidth', 1.5); % Second line in green
plot(x, mp_eps(3,:), 'b', 'LineWidth', 1.5); % Third line in blue
hold off; % Release the hold on the plot
saveas(h,"eg_supply_mp_3_0015-",'png')




