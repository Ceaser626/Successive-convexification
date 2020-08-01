function [x,u,sigma] =...
    Initialize_variable(k,K,m_wet,m_dry,r_I_init,v_I_init,v_I_fina,t_f_guess,g_I)

alpha1 = (K-k+1)/K;
alpha2 = (k-1)/K;

m = alpha1*m_wet+alpha2*m_dry;
r_I = alpha1*r_I_init;
v_I = alpha1*v_I_init+alpha2*v_I_fina;
q_B_I = [1,0,0,0]';
omiga_B = [0,0,0]';

x = [m,r_I',v_I',q_B_I',omiga_B']';
u = -m*g_I;
sigma = t_f_guess;