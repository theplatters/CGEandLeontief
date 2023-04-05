function [outineq, Out, outineq2, OutDeriv]=Derivs_HT_new_reallocI(X,A, beta, Omega, theta, sigma, L, N, N_orig, M, bar_X, kappa) 
 
p = X(1:N-M);
y = X(N-M+1:2*(N-M));

% Omega: (3N+M(3*N+1)) times  (3N+M(3*N+1)) factors -- the last M(N+1) rows
% L: a M(N+1) times 1 vector of labor endowments 
% theta -- (3N+M(3*N+1)) times 1 vector
% y -- 3N+M(3*N+1) vector
% p -- 3N+M(3*N+1) vector 
% Omega_ii = 1 for factors only here
% Cobb-Douglas between reallocatable factor and non-reallocatable one

beta = beta(1:N-M); 
theta = theta(1:N-M);
A = A(1:N-M);

kappa_f = kappa(end-N_orig+1:end);

w = 1;

Omega_q = Omega(1:N-M, 1:N-M);
Omega_p = Omega(1:N-M, :);
 
C = w'*L + y(1)*p(1)-y(2)*p(2)+p(N_orig+1+2:2+2*N_orig)'*y(2+N_orig+1:2+2*N_orig) - p(2+2*N_orig+1:2+3*N_orig)'*y(2+2*N_orig+1:2+3*N_orig);

exp_1 = (diag(p.^theta.*A.^(theta-1)))*(Omega_q.*((ones(N-M, N-M)*diag(p)).^(repmat(-theta, 1, N-M))));

exp_2 = (diag(A.^(theta-1))*(Omega_p.*(repmat(horzcat(p',w'),N-M, 1).^(repmat((1-theta), 1, N)))*ones(N, 1))).^(1./(1-theta));

adj_c_2 = p(2*N_orig+1+2:3*N_orig+2)./((1 - (kappa_f./2).*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1).^2) - (y(2*N_orig+1+2:3*N_orig+2)./bar_X).*kappa_f.*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1));
adj_exp_2 = p(2)/((1 - (kappa(1)/2)*(y(2) - 1)^2) - (y(2))*kappa(1)*(y(2) - 1));
exp_2(1) = adj_exp_2;
%exp_2(2) = 1;
exp_2(N_orig+1+2:2*N_orig+2) = adj_c_2;

Out(1:N-M) = (p - exp_2);
%sum(abs(Out(1:N-M)))
Out(N-M+1:2*(N-M)) = y' -y'*exp_1- beta'*(diag(p)^(-1))*C;
adj_costs_out = (y(N_orig+1+2:2*N_orig+2))' - (y(2*N_orig+1+2:3*N_orig+2).*(1 - (kappa_f./2).*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1).^2))';
Out(N-M+1+2+2*N_orig:N-M+3*N_orig+2) = adj_costs_out;
adj_costs_out1 = y(1) - (y(2)*(1 - (kappa(1)/2)*(y(2) - 1)^2));
Out(N-M+2) = adj_costs_out1;
%sum(abs(Out(N-M+1:2*(N-M))))
Out = Out';
outineq = [];
outineq2 = [];
 
%DWDP = (Omega_l'.*(1./repmat(L, 1, N-M))).*repmat((y)', M, 1);
DCDP = y.*vertcat(1, -1, zeros(N_orig, 1), ones(N_orig, 1), -1*ones(N_orig, 1)); %checked

%DWDY = ((Omega_l'.*(1./repmat(L, 1, N-M)))).*repmat((p)', M, 1);

DCDY = p.*vertcat(1, -1, zeros(N_orig, 1), ones(N_orig, 1), -1*ones(N_orig, 1));%checked

DOut1DP =  eye(N-M) - diag(A.^(-1+theta).*exp_2.^(theta))*(Omega_q.*(repmat((p'), N-M, 1).^(-repmat(theta, 1, N-M))));
DOut1DP(N_orig+1+2:2*N_orig+2, :) = horzcat(zeros(N_orig, N_orig+2), eye(N_orig), - diag(1./((1 - (kappa_f./2).*(y(2*N_orig+2+1:3*N_orig+2)./bar_X - 1).^2) - (y(2*N_orig+1+2:3*N_orig+2)./bar_X).*kappa_f.*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1))));
DOut1DP(1, :) = [1, -1/((1 - (kappa(1)/2)*(y(2) - 1)^2) - (y(2))*kappa(1)*(y(2) - 1)), zeros(1, N-M-2)];
%DOut1DP(2, :) = [0, 1, zeros(1, N-M-2)];

DOut1DY =  zeros(N-M, N-M);
DOut1DY(N_orig+1+2:2*N_orig+2, :) = horzcat(zeros(N_orig, 2*N_orig+2), diag((((adj_c_2).^2)./p(2*N_orig+1+2:3*N_orig+2)).*(- 2*(1./bar_X).*kappa_f.*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1)-(y(2*N_orig+1+2:3*N_orig+2)./((bar_X).^2)).*kappa_f )));
DOut1DY(1, :) = [0, p(2)*(-2*kappa(1)*(y(2) - 1)-(y(2))*kappa(1))/((1 - (kappa(1)/2)*(y(2) - 1)^2) - (y(2))*kappa(1)*(y(2) - 1))^2, zeros(1, N-M-2)];
%DOut1DY(2, :) = zeros(1, N-M);

DOut2DP = -(diag(theta)*diag(y)*diag(p)^(-1)*exp_1)'+diag(sum(diag(y.*p.^theta.*A.^(theta-1).*(theta))*(Omega_q.*(ones(N-M, N-M)*diag(p)).^(repmat(-theta-1, 1, N-M))), 1))...
    + 1*diag(beta.*p.^(-1-1))*C - bsxfun(@times, beta.*(p.^(-1)), DCDP');

DOut2DP(2*N_orig+1+2:3*N_orig+2, :) = zeros(N_orig, 3*N_orig+2);
DOut2DP(2, :) = zeros(1, N-M);

DOut2DY = eye(N-M) - exp_1' - bsxfun(@times, beta.*(p.^(-1)), DCDY');
DOut2DY(2*N_orig+1+2:3*N_orig+2, :) = horzcat(zeros(N_orig, N_orig+2), eye(N_orig), -diag((1 - (kappa_f./2).*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1).^2) - (y(2*N_orig+1+2:3*N_orig+2)./bar_X).*kappa_f.*(y(2*N_orig+1+2:3*N_orig+2)./bar_X - 1)));
DOut2DY(2, :) = [1, -((1 - (kappa(1)/2)*(y(2) - 1)^2) - (y(2))*kappa(1)*(y(2) - 1)), zeros(1, N-M-2)];

OutDeriv = [DOut1DP DOut1DY; DOut2DP DOut2DY]';
 
end