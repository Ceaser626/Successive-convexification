function [A_bar,B_bar,C_bar,Sigma_bar,Z_bar]=...
    Compute_transition(A,B,Sigma,z,K,k)

%------Compute A_bar
tao = 1/(K-1);
t_initial = (k-1)*tao;
t_final = k*tao;
A_bar = expm(A*tao);

%------Define derative             
B_bar = @(x)(expm(A.*(t_final-x)))*B.*(((t_final-x))/tao);
C_bar = @(x)(expm(A.*(t_final-x)))*B.*((x-t_initial)/tao);
Sigma_bar = @(x)(expm(A.*(t_final-x)))*Sigma;
Z_bar = @(x)(expm(A.*(t_final-x)))*z;

%------Compute B_bar、C_bar、Sigma_bar、Z_bar
B_bar = double(integral(B_bar,t_initial,t_final,'ArrayValued',true));
C_bar = double(integral(C_bar,t_initial,t_final,'ArrayValued',true));
Sigma_bar = double(integral(Sigma_bar,t_initial,t_final,'ArrayValued',true));
Z_bar = double(integral(Z_bar,t_initial,t_final,'ArrayValued',true));

% %------Compute A_bar
% tao = 1/(K-1);
% A_bar = expm(A*tao);
% 
% %------Define derative             
% B_bar = @(x)(expm(A.*(tao-x)))*B.*(1-x/tao);
% C_bar = @(x)(expm(A.*(tao-x)))*B.*(x/tao);
% Sigma_bar = @(x)(expm(A.*(tao-x)))*Sigma;
% Z_bar = @(x)(expm(A.*(tao-x)))*z;
% 
% %------Compute B_bar、C_bar、Sigma_bar、Z_bar
% B_bar = double(integral(B_bar,0,tao,'ArrayValued',true));
% C_bar = double(integral(C_bar,0,tao,'ArrayValued',true));
% Sigma_bar = double(integral(Sigma_bar,0,tao,'ArrayValued',true));
% Z_bar = double(integral(Z_bar,0,tao,'ArrayValued',true));



