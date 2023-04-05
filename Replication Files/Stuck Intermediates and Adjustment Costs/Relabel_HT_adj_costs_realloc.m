function [ NewOmega, Elasticities, FinalDemand] = Relabel_HT_adj_costs_realloc(Omega, alpha,epsilon,theta,beta, N, M,  sigma)
% N -- number of sectors
% M -- number of inputs

%Create a new Omega with first N rows: original industries, next N rows:
%corresponding intermediate inputs industries, next N rows: corresponding
%value added industries, last 2 rows: labour and capital.
%Then each industry has only one elasticity of substitution parameter put into Elasticities
%Finally, adjust final demand shares vector to have the correct size

%FShares N times M array of net factor shares 
%S_r N times M array: net share of a factor that cannot be reallocated

% Add M*N dummy non-reallocatable factors. S_r matters  
% Int_shares -- N^2 times 1 vector of shares of specific intermediate goods 
% Suppose this is the same across int sectors and final good sectors
% first intermediates for sector 1(all, in order) then, intermediates in
% sector 2, etc.

%%
NewOmega = zeros(3*N+3, 3*N+3);

for i=1:N
    NewOmega(i+2,N+i+2) = 1-alpha(i, 1); % original industries use their corresponding intermediate inputs...
    NewOmega(i+2,3*N+3) = alpha(i, 1); %... and their corresponding value added
end


NewOmega(N+1+2:2*N+2, 2*N+1+2:3*N+2) = eye(N);
NewOmega(2*N+1+2:3*N+2, 1+2:N+2) = Omega;
%%
Elasticities=[1; sigma; epsilon;ones(N, 1); theta; ones(M, 1)]; % Elasticities of substitution for transformed model
% epsilon is elasticity of substitution between VA and intermediates 
% N times 1
% theta is elasticity of substitution between intermediate good bundles 
% N times 1
% eta is elasticity of substitution between labour and capital -- no
% capital here so these are N times 1 ones 
% then, next N^2 industries are cd, so  -- ones as well as next N*(N+1)
% 


%%
FinalDemand=[zeros(2, 1); beta; zeros(M+2*N,1)]'; % no final demand for intermediate industries
NewOmega(2,:) = FinalDemand;
NewOmega(1, :) = [0; 1; zeros(3*N+M, 1)]';
FinalDemand=[1; zeros(3*N+M+1, 1)]; 
end

