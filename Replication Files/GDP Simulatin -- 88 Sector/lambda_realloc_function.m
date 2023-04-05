function lambda =lambda_realloc_function(TFP,init, beta, Omega, alpha, epsilon, theta, sigma,L) % no reallocation of labor
delta = 0.00001;
N = length(TFP);
Id = eye(N);
lambda = (beta'*inv(eye(N)-diag(1-alpha)*Omega))';
parfor i = 1:N
    x = TFP + Id(:,i)*delta;  
    lambda(i) = (GDP_realloc_function(x,init, beta, Omega, alpha, epsilon, theta, sigma,L)-GDP_realloc_function(TFP,init, beta, Omega, alpha, epsilon, theta, sigma,L))/delta;
end
lambda = lambda';

end 