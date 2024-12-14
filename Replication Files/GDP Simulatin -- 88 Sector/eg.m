function [GDP,mean_prices] = eg(shock,epsilon,theta,sigma,beta, Omega,alpha, L,N)
       optionsf = optimoptions('fmincon', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Algorithm','interior-point', 'TolCon', 1*10^(-10), 'TolFun', 10^(-10));
    
    init = [exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(shock));(beta'*inv(eye(N)-diag(1-alpha)*Omega))'./exp(-inv(eye(N)-diag(1-alpha)*Omega)*log(shock))];  % judicious choice of starting values

    [Soln,~,exitfl] = fmincon(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs(X, shock, beta, Omega, alpha, epsilon, theta, sigma,L),optionsf);

    if exitfl == 1 || exitfl == 2% solver no error (fmincon)
         GDP = ((shock.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(Soln(N+1:2*N).^(1/epsilon)).*(1./L).^(1/epsilon))' * L;
         mean_prices = sum(Soln(N+1:end) .* Soln(1:N)) / sum(Soln(N+1:end));
    else
        GDP = 1.0
        mean_prices = 1.0
    end



end