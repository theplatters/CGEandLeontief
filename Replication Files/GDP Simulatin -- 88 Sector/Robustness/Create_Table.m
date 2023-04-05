% loads simulation results for robustness and computes relevant statistics
% and writes them to a text file. 

theta_grid = [.001;0.2;0.99];
sigma_grid = [.8;0.9;0.99];
epsilon_grid = [.5;0.4;.6;0.99];

step_size = 4; %number of years for size of shocks.

fileID = fopen('Results.txt','w');

formatSpec = 'Full Reallocation\n';
fprintf(fileID, formatSpec);

for step_size = 1:3:4
formatSpec = 'Step size is %g\n';
fprintf(fileID, formatSpec, step_size);
    for count1 = 1:length(epsilon_grid)
        for count2 = 1:length(theta_grid)
            for count3 = 1:length(sigma_grid)
                
                epsilon = epsilon_grid(count1); % To match volatilities
                theta = theta_grid(count2); % Elasticity of substitution between intermediates
                sigma = sigma_grid(count3); % Elasticity of substitution between in consumption
                %cd(sprintf('%s%g%s', 'C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\',step_size,' year\Full Reallocation'))
                load(sprintf('%s%g%s%s%s%g%s%g%s%g%s','C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\',step_size,' year\Full Reallocation\', 'GDP','sigma',sigma,'epsilon', epsilon,'theta', theta,'.mat'))
                
                formatSpec = '(%g, %g, %g) * %g * %g * %g * %g * %g \n';
                Domar_vol = domar_weights'*std(bsxfun(@minus, log(lambda_simul(:,correct))',log(domar_weights)'))';
                %Domar_vol = domar_weights'*(std((log(lambda_simul(:,correct)'))))';
                fprintf(fileID, formatSpec, sigma, epsilon, theta, mean(log(GDP(correct))), std(log(GDP(correct))), skewness(log(GDP(correct))), kurtosis(log(GDP(correct)))-3, Domar_vol);

                
            end
        end
    end
end

formatSpec = 'No Reallocation\n';
fprintf(fileID, formatSpec);

for step_size = 1:3:4
formatSpec = 'Step size is %g\n';
fprintf(fileID, formatSpec, step_size);
    for count1 = 1:length(epsilon_grid)
        for count2 = 1:length(theta_grid)
            for count3 = 1:length(sigma_grid)
                
                epsilon = epsilon_grid(count1); % To match volatilities
                theta = theta_grid(count2); % Elasticity of substitution between intermediates
                sigma = sigma_grid(count3); % Elasticity of substitution between in consumption
                load(sprintf('%s%g%s%s%s%g%s%g%s%g%s','C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\',step_size,' year\No Reallocation\', 'GDP','sigma',sigma,'epsilon', epsilon,'theta', theta,'.mat'))
                
                formatSpec = '(%g, %g, %g) * %g * %g * %g * %g * %g \n';
                Domar_vol = domar_weights'*std(bsxfun(@minus, log(lambda_simul(:,correct))',log(domar_weights)'))';
                fprintf(fileID, formatSpec, sigma, epsilon, theta, mean(log(GDP(correct))), std(log(GDP(correct))), skewness(log(GDP(correct))), kurtosis(log(GDP(correct)))-3, Domar_vol);
               
                
            end
        end
    end
end


%%
% record domar weightrs over time
count  =1 ;
for year = 1960:2005
%year = 1960;
temp=80:1:88;
temp=[8;60;62;temp'];
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



formatSpec = 'Volatility of Domar weights for step size 1 is %g in the data\n';
Data_domar_vol = domar_weights'*(std(diff(log(lambda'),1)))';
fprintf(fileID, formatSpec, Data_domar_vol);

formatSpec = 'Volatility of Domar weights for step size 4 is %g in the data\n';
Data_domar_vol = domar_weights'*std(diff(log(lambda(:,1:4:end)'),1))';
fprintf(fileID, formatSpec, Data_domar_vol);

fclose(fileID);
%% Generate Histograms

figure
clear
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\1 year\No Reallocation')
load('GDPsigma0.9epsilon0.5theta0.001.mat', 'GDP','correct')
GDP_annual = GDP(correct); 
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\4 year\No Reallocation')
load('GDPsigma0.9epsilon0.5theta0.001.mat', 'GDP','correct')
GDP_quad = GDP(correct); 
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\4 year\Full Reallocation')
load('GDPsigma0.99epsilon0.99theta0.99.mat', 'GDP','correct')
GDP_quad_CD = GDP(correct); 
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\1 year\Full Reallocation')
load('GDPsigma0.99epsilon0.99theta0.99.mat', 'GDP','correct')
GDP_ann_CD = GDP(correct); 
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\1 year\No Reallocation')
load('GDPsigma0.7epsilon0.3theta0.001.mat', 'GDP','correct')
GDP_low_annual = GDP(correct); 
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\4 year\No Reallocation')
load('GDPsigma0.7epsilon0.3theta0.001.mat', 'GDP','correct')
GDP_low_quad = GDP(correct); 

close all

dist_ann = fitdist(GDP_annual,'kernel');
dist_low = fitdist(GDP_low_annual,'kernel');
dist_quad = fitdist(GDP_quad,'kernel');
dist_ann_CD= fitdist(GDP_ann_CD,'kernel');
dist_quad_CD= fitdist(GDP_quad_CD,'Kernel','kernel','epanechnikov');
dist_quad_low = fitdist(GDP_low_quad,'kernel');
xgrid = linspace(0.9,1.1,100);
hold on
plot(xgrid, pdf(dist_ann,xgrid),'-','LineWidth',3,'color','b')
plot(xgrid, pdf(dist_low,xgrid),'-*','LineWidth',3,'color','black')
plot(xgrid, pdf(dist_ann_CD,xgrid),'--','LineWidth',3,'color','r')
xlim([0.96,1.04])
ylabel('pdf')
legend('Benchmark Annual','Low Elasticity' ,'Cobb-Douglas', 'Location', 'NorthWest')
hold off
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\')
matlab2tikz('Histogram_annual.tex')

figure;
xgrid = linspace(0.9,1.1,100);
hold on
plot(xgrid, pdf(dist_quad,xgrid),'-','LineWidth',3,'color','b')
plot(xgrid, pdf(dist_quad_low,xgrid),'-*','LineWidth',3,'color','black')
plot(xgrid, pdf(dist_quad_CD,xgrid),'--','LineWidth',3,'color','r')
xlim([0.9,1.08])
ylabel('pdf')
legend('Benchmark Quadrennial','Low Elasticity', 'Cobb-Douglas', 'Location', 'NorthWest')
hold off
cd('C:\Users\Baqaee\Dropbox\Work\Hultens Theorem\Calibration\GDP Simulatin -- 88 Sector\Robustness\')
matlab2tikz('Histogram_quad.tex')

