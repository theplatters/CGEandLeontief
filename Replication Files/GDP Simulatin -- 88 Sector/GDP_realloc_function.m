function prices=GDP_realloc_function(TFP,init, beta, Omega, alpha, epsilon, theta, sigma,L) % no reallocation of labor

[Soln,~,exitfl] = knitromatlab(@(X) trivial(X),init,[],[],[],[],[],[], @(X)Simulation_Derivs_realloc(X,  TFP, beta, Omega, alpha, epsilon, theta, sigma,L),[],[],'Knitro_options.opt');
GDP = NaN;
N = length(TFP);
if exitfl == 0
      %GDP = (beta'*(Soln.^(1-sigma)))^(1/(sigma-1));
      prices = Soln;
            
end
            

end 