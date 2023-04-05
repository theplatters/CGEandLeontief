% loads simulation results for robustness with adjustment costs and computes relevant statistics
% and writes them to a text file. 

clear;
files = dir(fullfile('*.mat'));
N = length(files);

fileID = fopen('Results.txt','w');
fprintf(fileID, '$(\\sigma, \\theta,\\epsilon, \\kappa, t)$, $E(\\log(GDP))$, $std(\\log(GDP))$, $skewness(\\log(GDP))$, $kurtosis(\\log(GDP))-3$, std_\\, mean(lost_resources/GDP)) \n');
fprintf(fileID, '\n');
formatSpec = '(%g, %g, %g, %g, %g) * %g * %g * %g * %g * %g  * %g\n';
for i = 1:N
    load(files(i).name);    
    if length(L) >1
        fprintf(fileID, 'No reallocation \n');
        fprintf(fileID, formatSpec, sigma, epsilon(1), theta(1), kappa(5), step_size, mean(log(GDP(correct))), std(log(GDP(correct))), skewness(log(GDP(correct))), kurtosis(log(GDP(correct)))-3, domar_vol, mean(lost_res(correct)./GDP(correct)));
    else
        fprintf(fileID, 'Full reallocation \n');
        fprintf(fileID, formatSpec, sigma, epsilon(1), theta(1), kappa(5), step_size, mean(log(GDP(correct))), std(log(GDP(correct))), skewness(log(GDP(correct))), kurtosis(log(GDP(correct)))-3, domar_vol, mean(lost_res(correct)./GDP(correct)));
    end
end

fclose(fileID);