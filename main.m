%% main.m
% AE4320 Assignment 2026
% Aerodynamic Model Identification using the Two-Step Approach
% Tycho and Eltjo Assignment

clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Load one flight-data file
% -------------------------------------------------------------------------

defaultFile = fullfile('simData', 'da3211_1.mat');

if exist(defaultFile, 'file')
    dataFile = defaultFile;
else
    [file, path] = uigetfile('*.mat', 'Select one flight-data .mat file');

    if isequal(file, 0)
        error('No .mat file selected.');
    end

    dataFile = fullfile(path, file);
end

data = load(dataFile);

fprintf('Loaded file:\n%s\n\n', dataFile);

%% ------------------------------------------------------------------------
% 2. Extract variables
% -------------------------------------------------------------------------

t     = data.t(:);

alpha = data.alpha(:);
beta  = data.beta(:);

Ax = data.Ax(:);
Ay = data.Ay(:);
Az = data.Az(:);

p = data.p(:);
q = data.q(:);
r = data.r(:);

phi   = data.phi(:);
theta = data.theta(:);
psi   = data.psi(:);

vtas = data.vtas(:);

da = data.da(:);
de = data.de(:);
dr = data.dr(:);

Tc1 = data.Tc1(:);
Tc2 = data.Tc2(:);
Tc  = 0.5 * (Tc1 + Tc2);

u_n = data.u_n(:);
v_n = data.v_n(:);
w_n = data.w_n(:);

N  = length(t);
dt = mean(diff(t));

fprintf('Number of samples: %d\n', N);
fprintf('Average timestep: %.6f s\n\n', dt);

%% ------------------------------------------------------------------------
% 3. Compute body-axis airspeed components
% -------------------------------------------------------------------------

u_body = vtas .* cos(alpha) .* cos(beta);
v_body = vtas .* sin(beta);
w_body = vtas .* sin(alpha) .* cos(beta);

%% ------------------------------------------------------------------------
% 4. Generate approximate position from navigation-frame velocity
% -------------------------------------------------------------------------

x_E = zeros(N,1);
y_E = zeros(N,1);
z_E = zeros(N,1);

for k = 2:N
    dt_k = t(k) - t(k-1);

    x_E(k) = x_E(k-1) + u_n(k-1) * dt_k;
    y_E(k) = y_E(k-1) + v_n(k-1) * dt_k;
    z_E(k) = z_E(k-1) + w_n(k-1) * dt_k;
end

%% ------------------------------------------------------------------------
% 5. Plot angular rates
% -------------------------------------------------------------------------

figure;
plot(t, p, 'LineWidth', 1.2); hold on;
plot(t, q, 'LineWidth', 1.2);
plot(t, r, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Angular rate [rad/s]');
legend('p', 'q', 'r', 'Location', 'best');
title('Body angular rates');

%% ------------------------------------------------------------------------
% 6. Plot attitude angles
% -------------------------------------------------------------------------

figure;
plot(t, phi, 'LineWidth', 1.2); hold on;
plot(t, theta, 'LineWidth', 1.2);
plot(t, psi, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Attitude angle [rad]');
legend('\phi', '\theta', '\psi', 'Location', 'best');
title('Aircraft attitude angles');

%% ------------------------------------------------------------------------
% 7. Plot accelerations
% -------------------------------------------------------------------------

figure;
plot(t, Ax, 'LineWidth', 1.2); hold on;
plot(t, Ay, 'LineWidth', 1.2);
plot(t, Az, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
legend('A_x', 'A_y', 'A_z', 'Location', 'best');
title('Body accelerations');

%% ------------------------------------------------------------------------
% 8. Plot alpha and beta
% -------------------------------------------------------------------------

figure;
plot(t, alpha, 'LineWidth', 1.2); hold on;
plot(t, beta, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Angle [rad]');
legend('\alpha', '\beta', 'Location', 'best');
title('Angle of attack and sideslip angle');

%% ------------------------------------------------------------------------
% 9. Plot true airspeed
% -------------------------------------------------------------------------

figure;
plot(t, vtas, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('V_{TAS} [m/s]');
legend('V_{TAS}', 'Location', 'best');
title('True airspeed');

%% ------------------------------------------------------------------------
% 10. Plot throttle
% -------------------------------------------------------------------------

figure;
plot(t, Tc1, 'LineWidth', 1.2); hold on;
plot(t, Tc2, '--', 'LineWidth', 1.2);
plot(t, Tc, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Throttle coefficient [-]');
legend('T_{c1}', 'T_{c2}', 'Average T_c', 'Location', 'best');
title('Throttle inputs');

%% ------------------------------------------------------------------------
% 11. Plot control surface deflections
% -------------------------------------------------------------------------

figure;
plot(t, da, 'LineWidth', 1.2); hold on;
plot(t, de, 'LineWidth', 1.2);
plot(t, dr, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Deflection [rad]');
legend('\delta_a', '\delta_e', '\delta_r', 'Location', 'best');
title('Control surface deflections');

%% ------------------------------------------------------------------------
% 12. Plot body-axis airspeed components
% -------------------------------------------------------------------------

figure;
plot(t, u_body, 'LineWidth', 1.2); hold on;
plot(t, v_body, 'LineWidth', 1.2);
plot(t, w_body, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Body-axis airspeed [m/s]');
legend('u', 'v', 'w', 'Location', 'best');
title('Body-axis airspeed components');

%% ------------------------------------------------------------------------
% 13. Plot generated trajectory
% -------------------------------------------------------------------------

figure;
plot3(x_E, y_E, z_E, 'LineWidth', 1.2);
grid on;
xlabel('North position x_E [m]');
ylabel('East position y_E [m]');
zlabel('Down position z_E [m]');
title('Generated aircraft trajectory');

%% ------------------------------------------------------------------------
% 14. Store processed data in one structure
% -------------------------------------------------------------------------

proc = struct();

proc.dataFile = dataFile;

proc.t  = t;
proc.dt = dt;
proc.N  = N;

proc.alpha = alpha;
proc.beta  = beta;

proc.Ax = Ax;
proc.Ay = Ay;
proc.Az = Az;

proc.p = p;
proc.q = q;
proc.r = r;

proc.phi   = phi;
proc.theta = theta;
proc.psi   = psi;

proc.vtas = vtas;

proc.u_body = u_body;
proc.v_body = v_body;
proc.w_body = w_body;

proc.da = da;
proc.de = de;
proc.dr = dr;

proc.Tc1 = Tc1;
proc.Tc2 = Tc2;
proc.Tc  = Tc;

proc.u_n = u_n;
proc.v_n = v_n;
proc.w_n = w_n;

proc.x_E = x_E;
proc.y_E = y_E;
proc.z_E = z_E;

fprintf('main.m finished successfully.\n');