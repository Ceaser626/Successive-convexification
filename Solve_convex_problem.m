function [X,U,sigma,V,Delta,Delta_sigma]=...
    Solve_convex_problem(X_old,U_old,sigma_old,A_Bar,B_Bar,C_Bar,Sigma_Bar,...
    Z_Bar,omiga_v,omiga_Delta,omiga_Delta_sigma,x_init,x_fina,K,m_dry,...
    theta_max,omiga_max,T_min,T_max,gamma_gs,delta_max)

cvx_begin
%------Define variables
variable X(14,K);
variable U(3,K);   
variable time;
variable V(14*(K-1),1);
variable Delta(1,K); 
variable Delta_sigma;

%------Cost function
minimize (time+omiga_v*norm(V,1)+omiga_Delta*norm(Delta,2)+omiga_Delta_sigma*Delta_sigma);

subject to
%------Boundary constraints
       X(1:7,1) == x_init(1:7,1);
       X(12:14,1) == x_init(12:14,1);
       X(2:14,K) == x_fina(2:14,1);

%------Dynamics constraints
for k = 1:(K-1)
    X(:,(k+1)) == A_Bar(:,(14*k-13):14*k)*X(:,k)+B_Bar(:,(3*k-2):3*k)*U(:,k)+C_Bar(:,(3*k-2):3*k)*U(:,(k+1))+Sigma_Bar(:,k)*time+Z_Bar(:,k)+V((14*k-13):14*k,1);
end

%------State constraints
for k = 1:K
    X(1,k) >= m_dry;
    tan(gamma_gs)*norm(X(3:4,k),2) <= X(2,k);
    norm(X(10:11,k),2) <= sqrt((1-cos(theta_max))/2);
    norm(X(12:14,k),2) <= omiga_max;
end

%------Control constraints
for k = 1:K
    T_min <= U_old(:,k)'*U(:,k)/norm(U_old(:,k),2);
    norm(U(:,k),2) <= T_max;
    cos(delta_max)*norm(U(:,k),2) <= U(1,k);
end

%------Trust regions
       Delta_X = X-X_old;
       Delta_U = U-U_old;
       Delta_time = time-sigma_old;
for k = 1:K
    Delta_X(:,k)'*Delta_X(:,k)+Delta_U(:,k)'*Delta_U(:,k) <= Delta(:,k); 
end
    norm(Delta_time,1) <= Delta_sigma;      

cvx_end

X = double(X);
U = double(U);
sigma = double(time); 
V = double(V);
Delta = double(Delta);
Delta_sigma = double(Delta_sigma);









