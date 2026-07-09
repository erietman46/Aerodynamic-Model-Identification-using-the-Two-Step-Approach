clear; clc; close all;

%LOAD EVERY VARIABLE
fileName = 'C:\Users\eltjo\OneDrive - Delft University of Technology\TU Delft\Specialization\Aerodynamic-Model-Identification-using-the-Two-Step-Approach\simData\da3211_1.mat';
variables = load(fileName);

disp(fieldnames(variables));

t = variables.t(:);
u_n = variables.u_n(:);
v_n = variables.v_n(:);
w_n = variables.w_n(:);
phi = variables.phi(:);
theta = variables.theta(:);
psi = variables.psi(:);
Ax = variables.Ax(:);
Ay = variables.Ay(:);
Az = variables.Az(:);
p = variables.p(:);
q = variables.q(:);
r = variables.r(:);
vtas = variables.vtas(:);
alpha = variables.alpha(:);
beta  = variables.beta(:);

%WIND KNOWN VARIABLES
Wx = 6;
Wy = 2;
Wz = 10;

%NUMERICALLY INTEGRATE SUCH THAT THE EXACT KNOWN COORDINATES ARE DETERMINED
%WHICH ARE LATER USED FOR CREATING NOISY GPS MEASUREMENTS FOR THE KALMAN
%FILTER
xE = cumtrapz(t, u_n + Wx);
yE = cumtrapz(t, v_n + Wy);
zE = cumtrapz(t, w_n + Wz);

figure;
plot3(xE, yE, zE);
grid on;
xlabel('North position x_E [m]');
ylabel('East position y_E [m]');
zlabel('Altitude-like coordinate z_E [m]');
title('Generated aircraft trajectory for da3211\_1');

%CREATING NOISY IMU MEASUREMENTS
lambda_Ax = 0.02; %bias term
lambda_Ay = 0.02; %bias term
lambda_Az = 0.03; %bias term

lambda_p = deg2rad(0.005); %bias term
lambda_q = deg2rad(0.005); %bias term
lambda_r = deg2rad(0.002); %bias term

sigma_p = deg2rad(0.005); 
sigma_q = deg2rad(0.005);
sigma_r = deg2rad(0.002);

sigma_Ax = 0.02;
sigma_Ay = 0.02;
sigma_Az = 0.03;

Ax_m = Ax + lambda_Ax + sigma_Ax*randn(size(t));
Ay_m = Ay + lambda_Ay + sigma_Ay*randn(size(t));
Az_m = Az + lambda_Az + sigma_Az*randn(size(t));

p_m = p + lambda_p + sigma_p*randn(size(t));
q_m = q + lambda_q + sigma_q*randn(size(t));
r_m = r + lambda_r + sigma_r*randn(size(t));

%CREATING THE NOISY INPUT VECTOR u_m
u_IMU = [Ax_m, Ay_m, Az_m, p_m, q_m, r_m];

%CREATING NOISY GPS MEASUREMENTS
sigma_pos = 1.0;              % m
sigma_vel = 0.01;             % m/s
sigma_att = deg2rad(0.04);    % rad

x_GPS = xE + sigma_pos*randn(size(t));
y_GPS = yE + sigma_pos*randn(size(t));
z_GPS = zE + sigma_pos*randn(size(t));

u_GPS = u_n + Wx + sigma_vel*randn(size(t));
v_GPS = v_n + Wy + sigma_vel*randn(size(t));
w_GPS = w_n + Wz + sigma_vel*randn(size(t));

phi_GPS   = phi   + sigma_att*randn(size(t));
theta_GPS = theta + sigma_att*randn(size(t));
psi_GPS   = psi   + sigma_att*randn(size(t));

%CREATING NOISY AIRDATA MEASUREMENTS
sigma_V = 0.2;
sigma_alpha = deg2rad(0.1);
sigma_beta  = deg2rad(0.25);

V_m     = vtas  + sigma_V*randn(size(t));
alpha_m = alpha + sigma_alpha*randn(size(t));
beta_m  = beta  + sigma_beta*randn(size(t));

%CREATING THE NOISY OUTPUT VECTOR z_m
z_m = [x_GPS, y_GPS, z_GPS, ...
       u_GPS, v_GPS, w_GPS, ...
       phi_GPS, theta_GPS, psi_GPS, ...
       V_m, alpha_m, beta_m];







