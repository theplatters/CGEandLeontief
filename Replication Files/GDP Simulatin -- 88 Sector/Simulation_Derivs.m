function [outineq, Out, outineq2, OutDeriv]=Simulation_Derivs(X, A, beta, Omega, alpha, epsilon, theta, sigma,L) % no reallocation of labor
N = length(alpha);
p = X(1:N);
y = X(N+1:2*N);

q = (Omega*(p.^(1-theta))).^(1/(1-theta));
w = p.*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon)).*(1./L).^(1/epsilon);
C = w'*L;

Out(1:N) = p - (diag(A)^(epsilon-1)*(alpha.*w.^(1-epsilon)+(1-alpha).*q.^(1-epsilon))).^(1/(1-epsilon));
Out(N+1:2*N) = y' -y'*diag(p)^epsilon*diag(A)^(epsilon-1)*diag(q)^(theta-epsilon)*diag(1-alpha)* Omega * diag(p)^(-theta) - beta'*diag(p)^(-sigma)*C;
Out = Out';
outineq = [];
outineq2 = [];

DQDP = bsxfun(@times, (q.^theta), (p.^(-theta))').*Omega; %
DWDP = diag((A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon)).*(1./L).^(1/epsilon)); %checked
DWDY = (1/epsilon)*diag(p.*(A.^((epsilon-1)/epsilon)).*(alpha.^(1/epsilon)).*(y.^(1/epsilon-1)).*(L).^(-1/epsilon)); %checked
DCDP = DWDP'*L; %checked
DCDY = DWDY'*L;%checked
DOut1DP = eye(N) - diag(diag(A)^(-1)*((alpha.*(w.^(1-epsilon))+(1-alpha).*(q.^(1-epsilon)))).^(epsilon/(1-epsilon)))*...
    (diag(alpha)*diag(w.^(-epsilon))*DWDP+diag(1-alpha)*diag(q.^(-epsilon))*DQDP);

DOut1DY =  -diag(diag(A.^(-1))*((alpha.*(w.^(1-epsilon))+(1-alpha).*(q.^(1-epsilon)))).^(epsilon/(1-epsilon)))*...
    (diag(alpha)*diag(w.^(-epsilon))*DWDY);

DOut2DP = -(epsilon * diag(p.^(-theta))*Omega'*diag((p.^(epsilon-1)).*(y).*(q.^(theta-epsilon)).*(1-alpha).*(A.^(epsilon-1)))...
    +(theta-epsilon)*diag(p.^(-theta))*Omega'*diag((p.^(epsilon)).*(y).*(q.^(theta-epsilon-1)).*(1-alpha).*(A.^(epsilon-1)))*DQDP ...
    -sigma*diag(beta.*p.^(-sigma-1))*C+ bsxfun(@times, beta.*(p.^(-sigma)), DCDP')...
    - theta* diag(p.^(-theta-1)).*diag(Omega'*diag((p.^(epsilon)).*(q.^(theta-epsilon)).*(1-alpha).*(A.^(epsilon-1)))*y));

DOut2DY = eye(N) - (diag(p)^epsilon*diag(A)^(epsilon-1)*diag(q)^(theta-epsilon)*diag(1-alpha)*Omega*diag(p)^(-theta))' - bsxfun(@times, beta.*(p.^(-sigma)), DCDY');
OutDeriv = [DOut1DP DOut1DY; DOut2DP DOut2DY]';

end
