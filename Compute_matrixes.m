function [A_bar,B_bar,C_bar,Sigma_bar,Z_bar]=...
    Compute_matrixes(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B,K,k)

%------Compute A、B、Sigma and z
%---Call function 'Compute_formula'
[A,B,Sigma,z]=...
    Compute_formula(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B);

%------Compute A_bar、B_bar、C_bar、Sigma_bar and Z_bar
%---Call function 'Compute_transition'       
[A_bar,B_bar,C_bar,Sigma_bar,Z_bar]=...
    Compute_transition(A,B,Sigma,z,K,k);

