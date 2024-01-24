close all; clear; clc;

% helicopter basic parameters
R = 4.4102;
Omega = 43.6825;
b = 3;
e = 0.05*R;
Mb = 110.75/b;

% uniform blade assumption
S_beta = Mb*(R - e)/2;
I_beta = Mb*(R - e)^2/3;
K_beta = 0.;

% global parameters (used later to compute function and Jacobian matrix)
global h x_CM ell K_H A v_tip sigma C_Lalpha C_D f S_bar ell_bar C_m_f C_L_calpha S_c W rho

h = 2.; % m
x_CM = -0.24; % m
ell = 6.; % m
K_H = b/2*(e*S_beta*Omega^2 + K_beta); % N m/radian
A = pi*R^2; % m^2
v_tip = R*Omega; % m/s
sigma = 0.10195;
% C_Lalpha = 5.73; % typical coefficient (91% of 2*pi)
% C_Lalpha = 6.0447; % NACA 23012 at Mach 0-0.2 used in MBDyn
C_Lalpha = 7.0760; % NACA 23012 at Mach 0.5 used in MBDyn
C_D = 0.008;
f = .4; % m^2
S_bar = 2.; % m^2
ell_bar = 1.; % m
C_m_f = 0.02;
C_L_calpha = 5.73;
S_c = 1.2; % m^2
W = 9.8*1190; % N
rho = 1.1116; % kg/m^3

% initialize variables (initial guess)
T_D = W; % N
H_D = W/20; % N
R_f = W/20; % N
M_f = 0.; % Nm
P_c = 1000.; % N
a_1 = 0.; % radian
theta_0 = 8/180*pi; % radian
B_1 = 0.; % radian
gamma = 0.; % radian
% tau = -0.2:0.02:0.2; % 0.; % radian (inverse problem) fix
% V_infty = 50; % 0:5:100; % m/s (inverse problem) fix
tau = 0.; % radian (inverse problem) fix
V_infty = 0:5:75; % m/s (inverse problem) fix
mu = 0/v_tip;
lambda = 0.02;
alpha_D = 0.; % radian
u = lambda*v_tip; % m/s
v_1 = v_tip*sqrt(mu^2 + lambda^2); % m/s

vars_inverse(:,1) = [T_D;H_D;a_1;theta_0;B_1;gamma;R_f;M_f;P_c;mu;lambda;alpha_D;u;v_1];
% n_var = tau;
n_var = V_infty;
vars_inverse_loop = zeros(14,length(n_var));

for j = 1:length(n_var)

imax = 10;
tol = 1e-3;
i = 0;
err = 1;

while i<imax && err>tol

    i = i + 1;
    [F,J] = trim_jacobian_inverse(Omega, R, V_infty(j), tau, vars_inverse(:,i)); % DON'T forget to change the j position!!!
    vars_inverse(:,i+1) = vars_inverse(:,i) - J\F;

    err = norm(F);
end
vars_inverse_loop(:,j) = vars_inverse(:,end);
end

 
alpha_H = (tau+vars_inverse_loop(6,:));
a_1_H = vars_inverse_loop(3,:) - vars_inverse_loop(5,:);
T_H = vars_inverse_loop(1,:).*cos(a_1_H) - vars_inverse_loop(2,:).*sin(a_1_H);
H_H = vars_inverse_loop(1,:).*sin(a_1_H) + vars_inverse_loop(2,:).*cos(a_1_H);

% Plots

figure(1)
plot(n_var, vars_inverse_loop(1,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('T_D (N)')

figure(2)
plot(n_var, vars_inverse_loop(2,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('H_D (N)')

figure(3)
plot(n_var, vars_inverse_loop(3,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('a_1 (radian)')


figure(4)
plot(n_var, vars_inverse_loop(4,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\theta_0 (radian)')

figure(5)
plot(n_var, vars_inverse_loop(5,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('B_1 (radian)')

figure(6)
plot(n_var, vars_inverse_loop(6,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\gamma (radian)')

figure(7)
plot(n_var, vars_inverse_loop(7,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('R_f (N)')

figure(8)
plot(n_var, vars_inverse_loop(8,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('M_f (Nm)')

figure(9)
plot(n_var, vars_inverse_loop(9,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('P_c (N)')

figure(10)
plot(n_var, vars_inverse_loop(10,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\mu')

figure(11)
plot(n_var, vars_inverse_loop(11,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\lambda')

figure(12)
plot(n_var, vars_inverse_loop(12,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\alpha_D (radian)')

figure(13)
plot(n_var, vars_inverse_loop(13,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('u (m/s)')

figure(14)
plot(n_var, vars_inverse_loop(14,:),'-.ob')
xlabel('V_\infty (m/s)')
ylabel('v_1 (m/s)')

figure(15)
plot(n_var, alpha_H,'-.ob')
xlabel('V_\infty (m/s)')
ylabel('\alpha_H (radian)')

figure(16)
plot(n_var, a_1_H,'-.ob')
xlabel('V_\infty (m/s)')
ylabel('a_{1H} (radian)')

figure(17)
plot(n_var, T_H,'-.ob')
xlabel('V_\infty (m/s)')
ylabel('T_H (N)')

figure(18)
plot(n_var, H_H,'-.ob')
xlabel('V_\infty (m/s)')
ylabel('H_H (N)')

%%
function [F,J] = trim_jacobian_inverse(Omega, R, V_infty, tau, vars)

global h x_CM ell K_H A v_tip sigma C_Lalpha C_D f S_bar ell_bar C_m_f C_L_calpha S_c W rho

T_D = sym('T_D');
H_D = sym('H_D');
a_1 = sym('a_1');
theta_0 = sym('theta_0');
B_1 = sym('B_1');
gamma = sym('gamma');
R_f = sym('R_f');
M_f = sym('M_f');
P_c = sym('P_c');
mu = sym('mu');
lambda = sym('lambda');
alpha_D = sym('alpha_D');
u = sym('u');
v_1 = sym('v_1');

% Equation
eq1 = T_D*cos(a_1-B_1)- H_D*sin(a_1-B_1) - W*cos(gamma) - R_f*sin(gamma+tau) + P_c*cos(gamma+tau);
eq2 = T_D*sin(a_1-B_1) + H_D*cos(a_1-B_1) - W*sin(gamma) + R_f*cos(gamma+tau) + P_c*sin(gamma+tau);
eq3 = K_H*(a_1-B_1) + W*h*sin(gamma) - R_f*h*cos(gamma+tau) + M_f - W*x_CM*cos(gamma) - R_f*x_CM*sin(gamma+tau) - P_c*ell*cos(gamma+tau);
eq4 = - T_D + rho*sigma*A*v_tip^2*C_Lalpha*0.5*(theta_0*(1/3 + mu^2/2) - lambda/2 - mu*a_1/2);
eq5 = - H_D + rho*sigma*A*v_tip^2*(mu*C_D/4 + C_Lalpha*1/4*(mu*lambda*theta_0 - lambda*a_1/2));
eq6 = - R_f + 0.5*rho*V_infty^2*f; 
eq7 = - M_f + 0.5*rho*V_infty^2*C_m_f*S_bar*ell_bar;
eq8 = - P_c - 0.5*rho*V_infty^2*C_L_calpha*S_c*(gamma + tau);
eq9 = - a_1 + (1/(1-0.5*mu^2))*(2*mu*(4/3 * theta_0 - lambda));
eq10 = - mu + (V_infty*cos(alpha_D))/(Omega*R);
eq11 = - lambda + (V_infty*sin(alpha_D)+u)/(Omega*R);
eq12 = - u + T_D/(2*rho*A*abs(v_1));
eq13 = - v_1 + Omega*R*sqrt((mu^2 + lambda^2));
eq14 = - alpha_D + tau + gamma - a_1 + B_1;

% result = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14);
J = jacobian([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14],[T_D,H_D,a_1,theta_0,B_1,gamma,R_f,M_f,P_c,mu,lambda,alpha_D,u,v_1]);
F = [eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11;eq12;eq13;eq14];

J = subs(J,{T_D;H_D;a_1;theta_0;B_1;gamma;R_f;M_f;P_c;mu;lambda;alpha_D;u;v_1},{vars});
F = subs(F,{T_D;H_D;a_1;theta_0;B_1;gamma;R_f;M_f;P_c;mu;lambda;alpha_D;u;v_1},{vars});

J = double(J);
F = double(F);
end