function [outineq, Out, outineq2, OutDeriv]=Simulation_Derivs_realloc(X, A, beta, Omega, alpha, epsilon, theta, sigma,L) % full reallocation
N = length(alpha);
p = X(1:N);

q = (Omega*(p.^(1-theta))).^(1/(1-theta));
w = 1;
%C = sum(L)/Pc;

Out(1:N) = p - (diag(A)^(epsilon-1)*(alpha.*w.^(1-epsilon)+(1-alpha).* q.^(1-epsilon))).^(1/(1-epsilon));
Out = Out';
outineq = [];
outineq2 = [];

DQDP = bsxfun(@times, (q.^theta), (p.^(-theta))').*Omega; % 
%DWDP = diag((A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon)).*(1./L).^(1/epsilon)); %checked
%DWDY = (1/epsilon)*diag(p.*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon-1)).*(L).^(-1/epsilon)); %checked
DOut1DP = eye(N) - diag(diag(A)^(-1)*((alpha.*(w.^(1-epsilon))+(1-alpha).*(q.^(1-epsilon)))).^(epsilon/(1-epsilon)))*...
    (diag(1-alpha)*diag(q.^(-epsilon))*DQDP);

OutDeriv = [DOut1DP]';

end

