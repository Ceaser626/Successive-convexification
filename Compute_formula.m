function [A,B,Sigma,z]=...
    Compute_formula(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B)

%------Define symbolic variables
syms m real
r_I = sym('r_I',[3 1]);
v_I = sym('v_I',[3 1]);
q_B_I = sym('q_B_I',[4 1]);
omiga_B = sym('omiga_B',[3 1]);
T_B = sym('T_B',[3 1]);

%------Define matrixes
xl = [m,r_I',v_I',q_B_I',omiga_B']';
ul = T_B;
C_B_I = [    1-2*(q_B_I(3)^2+q_B_I(4)^2)           2*(q_B_I(2)*q_B_I(3)+q_B_I(1)*q_B_I(4))   2*(q_B_I(2)*q_B_I(4)-q_B_I(1)*q_B_I(3));
         2*(q_B_I(2)*q_B_I(3)-q_B_I(1)*q_B_I(4))       1-2*(q_B_I(2)^2+q_B_I(4)^2)           2*(q_B_I(3)*q_B_I(4)+q_B_I(1)*q_B_I(2));
         2*(q_B_I(2)*q_B_I(4)+q_B_I(1)*q_B_I(3))   2*(q_B_I(3)*q_B_I(4)-q_B_I(1)*q_B_I(2))       1-2*(q_B_I(2)^2+q_B_I(3)^2)        ;];
C_I_B = C_B_I';
Omiga = [    0        -omiga_B(1)   -omiga_B(2)   -omiga_B(3);
         omiga_B(1)       0          omiga_B(3)   -omiga_B(2);
         omiga_B(2)   -omiga_B(3)       0          omiga_B(1);
         omiga_B(3)    omiga_B(2)   -omiga_B(1)        0     ;];
r_T_B_x = [    0      -r_T_B(3)  r_T_B(2);
            r_T_B(3)     0      -r_T_B(1);
           -r_T_B(2)   r_T_B(1)     0    ;];
omiga_B_x = [    0        -omiga_B(3)   omiga_B(2);
             omiga_B(3)       0        -omiga_B(1);
            -omiga_B(2)    omiga_B(1)       0     ;];
        
%------Define symbolic equation elements
m_dot = -alpha_m_dot*norm(T_B,2);
r_I_dot = v_I;
v_I_dot = 1/m*C_I_B*T_B+g_I;
q_B_I_dot = 1/2*Omiga*q_B_I;
omiga_B_dot = (J_B)\(r_T_B_x*T_B-omiga_B_x*J_B*omiga_B);

%------f(x,u)
f = [m_dot,r_I_dot',v_I_dot',q_B_I_dot',omiga_B_dot']';

%------Jacobian
jacobian_x = jacobian(f,xl);
jacobian_u = jacobian(f,ul);

%------Compute f、Jacobian
Symbol = [xl',ul']';
Number = [x',u']';
f = double(subs(f,Symbol,Number));
jacobian_x = double(subs(jacobian_x,Symbol,Number));
jacobian_u = double(subs(jacobian_u,Symbol,Number));

%------Compute A、B、Sigma、z
A = sigma*jacobian_x;
B = sigma*jacobian_u;
Sigma = f;
z = -A*x-B*u;





% function [A,B,Sigma,z]=...
%     Compute_formula(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B)
% 
% %------Define symbolic variables
% syms m real
% r_I = sym('r_I',[3 1]);
% v_I = sym('v_I',[3 1]);
% q_B_I = sym('q_B_I',[4 1]);
% omiga_B = sym('omiga_B',[3 1]);
% T_B = sym('T_B',[3 1]);
% 
% %------Define matrixes
% xl = [m,r_I',v_I',q_B_I',omiga_B']';
% ul = T_B;
% C_B_I = [    1-2*(q_B_I(3)^2+q_B_I(4)^2)           2*(q_B_I(2)*q_B_I(3)+q_B_I(1)*q_B_I(4))   2*(q_B_I(2)*q_B_I(4)-q_B_I(1)*q_B_I(3));
%          2*(q_B_I(2)*q_B_I(3)-q_B_I(1)*q_B_I(4))       1-2*(q_B_I(2)^2+q_B_I(4)^2)           2*(q_B_I(3)*q_B_I(4)+q_B_I(1)*q_B_I(2));
%          2*(q_B_I(2)*q_B_I(4)+q_B_I(1)*q_B_I(3))   2*(q_B_I(3)*q_B_I(4)-q_B_I(1)*q_B_I(2))       1-2*(q_B_I(2)^2+q_B_I(3)^2)        ;];
% C_I_B = C_B_I';
% Omiga = [    0        -omiga_B(1)   -omiga_B(2)   -omiga_B(3);
%          omiga_B(1)       0          omiga_B(3)   -omiga_B(2);
%          omiga_B(2)   -omiga_B(3)       0          omiga_B(1);
%          omiga_B(3)    omiga_B(2)   -omiga_B(1)        0     ;];
% r_T_B_x = [    0      -r_T_B(3)  r_T_B(2);
%             r_T_B(3)     0      -r_T_B(1);
%            -r_T_B(2)   r_T_B(1)     0    ;];
% omiga_B_x = [    0        -omiga_B(3)   omiga_B(2);
%              omiga_B(3)       0        -omiga_B(1);
%             -omiga_B(2)    omiga_B(1)       0     ;];
%         
% %------Define symbolic equation elements
% m_dot = -alpha_m_dot*norm(T_B,2);
% r_I_dot = v_I;
% v_I_dot = 1/m*C_I_B*T_B+g_I;
% q_B_I_dot = 1/2*Omiga*q_B_I;
% omiga_B_dot = (J_B)\(r_T_B_x*T_B-omiga_B_x*J_B*omiga_B);
% 
% %------f(x,u)
% f = [m_dot,r_I_dot',v_I_dot',q_B_I_dot',omiga_B_dot']';
% 
% %------Jacobian
% jacobian_x = jacobian(f,xl);
% jacobian_u = jacobian(f,ul);
% 
% %------Compute f、Jacobian
% Symbol = [xl',ul']';
% Number = [x',u']';
% f = double(subs(f,Symbol,Number));
% jacobian_x = double(subs(jacobian_x,Symbol,Number));
% jacobian_u = double(subs(jacobian_u,Symbol,Number));
% 
% %------Compute A、B、Sigma、z
% A = sigma*jacobian_x;
% B = sigma*jacobian_u;
% Sigma = f;
% z = -A*x-B*u;






