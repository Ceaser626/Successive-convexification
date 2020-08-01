%%Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time

%-----------------------------Initialization-------------------------------
clear all
clc
%------0.Set Parameters          
%---Model Parameters
g_I = [-1,0,0]';
m_wet = 2;
m_dry = 1;
T_max = 5;
T_min = 0.3;
I_sp = 10;                   %ungiven in paper
delta_max = 20/180*pi;
theta_max = 90/180*pi;
gamma_gs = 20/180*pi;
omiga_max = 60/180*pi;
J_B = 1e-2*eye(3);
r_T_B = -1e-2*[1,0,0]';
%---Algorithm Parameters
omiga_v = 1e5;
omiga_Delta = 1e-3;
omiga_Delta_sigma = 1e-1;
Delta_tol = 1e-3;
sigma_tol = 1e-3;           %ungiven in paper
v_tol = 1e-10;
N_iter_max = 25;
K = 50;
r_I_init = [4,4,0]';
r_I_fina = [0,0,0]';
% v_I_init = [0,-1,-2]';
v_I_init = [0,-0.75,0.5]';%ungiven in paper
v_I_fina = -1e-1*[1,0,0]';
q_B_I_init = [1,0,0,0]';    %ungiven in paper
q_B_I_fina = [1,0,0,0]';
omiga_B_init = [0,0,0]';
omiga_B_fina = [0,0,0]';    %ungiven in paper
t_f_guess = 3;              %ungiven in paper   
%---Compute Parameters
tao = 1/(K-1);
alpha_m_dot = 1/I_sp/abs(g_I(1))/10;
x_init = [m_wet,r_I_init',v_I_init',q_B_I_init',omiga_B_init']';
x_fina = [0,r_I_fina',v_I_fina',q_B_I_fina',omiga_B_fina']';
%---Define matrixes
A_Bar = zeros(14,14*(K-1));
B_Bar = zeros(14,3*(K-1));
C_Bar = zeros(14,3*(K-1));
Sigma_Bar = zeros(14,(K-1));
Z_Bar = zeros(14,(K-1));
X = zeros(14,K);
U = zeros(3,K);

%------1.First Iteration            
for k = 1:K
%---Call function 'Initialize_variable'
[x,u,sigma] =...
    Initialize_variable(k,K,m_wet,m_dry,r_I_init,v_I_init,v_I_fina,t_f_guess,g_I);
%---Store data
X(:,k) = x;
U(:,k) = u;
end

%------2.Second Iteration
for k = 1:(K-1)
%---Acquire x、u
x = X(:,k);
u = U(:,k);
%---Call function 'Compute_matrixes'
[A_bar,B_bar,C_bar,Sigma_bar,Z_bar]=...
    Compute_matrixes(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B,K,k);
%---Store data
A_Bar(:,(14*k-13):14*k) = A_bar;
B_Bar(:,(3*k-2):3*k) = B_bar;
C_Bar(:,(3*k-2):3*k) = C_bar;
Sigma_Bar(:,k) = Sigma_bar;
Z_Bar(:,k) = Z_bar;
end


%------------------------Convex Optimization Loop--------------------------

for iterate = 1:N_iter_max
    
%------Display Loop Number
display = ['Loop number ',num2str(iterate)];
disp(display)  

%------Solve Convexified Problem        
%---Call function 'Solve_convex_problem' 
disp('Solve convex problem')
[X,U,sigma,V,Delta,Delta_sigma]=...
    Solve_convex_problem(X,U,sigma,A_Bar,B_Bar,C_Bar,Sigma_Bar,...
    Z_Bar,omiga_v,omiga_Delta,omiga_Delta_sigma,x_init,x_fina,K,m_dry,...
    theta_max,omiga_max,T_min,T_max,gamma_gs,delta_max);
sigma_store(iterate,1) = sigma;

%------Judge Iteration Conditions
Diverse(1,iterate) = norm(Delta,2);
Diverse(2,iterate) = norm(V,1);
Diverse(3,iterate) = Delta_sigma;
if (norm(Delta,2) <= Delta_tol)&&(norm(V,1) <= v_tol)&&(Delta_sigma <=sigma_tol)
    break
end

%------Update A_Bar、B_Bar、C_Bar、Sigma_Bar、Z_Bar
for k = 1:(K-1)
%---Acquire x、u
x = X(:,k);
u = U(:,k);
%---Call function 'Compute_matrixes'
[A_bar,B_bar,C_bar,Sigma_bar,Z_bar]=...
    Compute_matrixes(x,u,sigma,alpha_m_dot,g_I,J_B,r_T_B,K,k);
%---Store data
A_Bar(:,(14*k-13):14*k) = A_bar;
B_Bar(:,(3*k-2):3*k) = B_bar;
C_Bar(:,(3*k-2):3*k) = C_bar;
Sigma_Bar(:,k) = Sigma_bar;
Z_Bar(:,k) = Z_bar;
end    
 
end


%--------------------------------Plot Figure-------------------------------

%------State plot
figure
%---Position--Time
display_position(1:3,:) = X(2:4,:);
subplot(3,3,1)
plot(0:tao:1,display_position)
grid on
legend('x','y','z');
xlabel('时间')
ylabel('位置')
title('Position--Time')
%---Velocity--Time
display_velocity(1:3,:) = X(5:7,:);
subplot(3,3,2)
plot(0:tao:1,display_velocity)
grid on
legend('x','y','z');
xlabel('时间')
ylabel('速度')
title('Velocity--Time')
%---Thrust--Time
display_netthrust = U;
subplot(3,3,3)
plot(0:tao:1,display_netthrust)
grid on
legend('x','y','z');
xlabel('时间')
ylabel('净推力')
title('Thrust--Time')
%---Mass--Time
mass = X(1,:);
subplot(3,3,4)
plot(0:tao:1,mass)
grid on
xlabel('时间')
ylabel('质量')
title('Mass--Time')
%---Allvelocity--Time
display_allvelocity = zeros(1,K);
for i = 1:K
display_allvelocity(i) = norm(X(5:7,i),2);
end
subplot(3,3,5)
plot(0:tao:1,display_allvelocity)
grid on
xlabel('时间')
ylabel('总速度')
title('Allvelocity--Time')
%---Allthrust--Time
display_allthrust = zeros(1,K);
for i = 1:K
display_allthrust(i) = norm(U(:,i),2);
end
subplot(3,3,6)
hold on
plot(0:tao:1,display_allthrust,'-')
plot(0:tao:1,display_allthrust,'.')
grid on
xlabel('时间')
ylabel('总净推力')
title('Allthrust--Time')
hold off
%---omiga--Time
display_omiga(1:3,:) = X(12:14,:);
subplot(3,3,7)
plot(0:tao:1,display_omiga)
grid on
legend('x','y','z');
xlabel('omiga_x')
ylabel('omiga_y')
zlabel('omiga_z')
title('omiga--Time')
%---Allomiga--Time
display_allomiga = zeros(1,K);
for i = 1:K
display_allomiga(i) = norm(X(12:14,i),2)*180/pi;
end
subplot(3,3,8)
plot(0:tao:1,display_allomiga)
grid on
xlabel('时间')
ylabel('总角速度')
title('Allomiga--Time')
%---Tilt--Time
display_tilt = zeros(1,K);
for i = 1:K
display_tilt(i) = acos(1-2*(norm(X(10:11,i),2))^2)/pi*180;
end
subplot(3,3,9)
plot(0:tao:1,display_tilt)
grid on
xlabel('时间')
ylabel('倾斜角')
title('Tilt--Time')

%------Trajectory plot
figure
%---North--Up
subplot(2,2,4)
plot(X(4, :), X(2, :))
grid on
xlabel('North-Position')  
ylabel('Up-Position')
title('North--Up')
%---East--Up
subplot(2,2,3)
plot(X(3, :), X(2, :))
grid on
xlabel('East-Position')  
ylabel('Up-Position')
title('East--Up')
%---3D Trajectory
subplot(2,2,2)
plot3(X(3,:),X(4,:),X(2,:))
grid on
xlabel('East-Position')  
ylabel('North-Position')
zlabel('Up-Position')
title('3D Trajectory')
view(-142.5,30);
% set(gca,'DataAspectRatio',[1 1 1])
%---East--North
subplot(2,2,1)
plot(X(3, :), X(4, :))
grid on
xlabel('East-Position')  
ylabel('North-Position')
title('East--North')

%------Iterate plot
figure
%---Error state & control
subplot(2,2,1)
plot(1:iterate,Diverse(1,:))
grid on
xlabel('迭代次数')
ylabel('状态&控制误差')
title('state & control error--iterate time')
%---Error visual control
subplot(2,2,2)
plot(1:iterate,Diverse(2,:))
grid on
xlabel('迭代次数')
ylabel('虚拟控制')
title('visual control--iterate time')
%---Error time
subplot(2,2,3)
plot(1:iterate,Diverse(3,:))
grid on
xlabel('迭代次数')
ylabel('估计时间')
title('flight time error--iterate time')
%---Flight time
subplot(2,2,4)
plot(1:iterate,sigma_store)
grid on
xlabel('迭代次数')
ylabel('飞行时间')
title('flight time--iterate time')


