function [outineq, Out]=Simulation(X,A, beta, Omega, alpha, epsilon, theta, sigma,L)
N = length(alpha);
p = X(1:N);
y = X(N+1:2*N);

q = (Omega*(p.^(1-theta))).^(1/(1-theta));
w = p.*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon)).*(1./L).^(1/epsilon);
C = w'*L;

Out(1:N) = p - (diag(A)^(epsilon-1)*(alpha.*w.^(1-epsilon)+(1-alpha).*q.^(1-epsilon))).^(1/(1-epsilon));
Out(N+1:2*N) = y' -y'*diag(p)^epsilon*diag(A)^(epsilon-1)*diag(q)^(theta-epsilon)*diag(1-alpha)*Omega*diag(p)^(-theta) - beta'*diag(p)^(-sigma)*C;
Out = Out';
outineq = 0;

end