clc;
clear;
close all;

% Simple 18-state Extended Kalman Filter for maneuver DA3211_1

rng(2026);
g = 9.803;

%% Load the data made in Part 1

load('part1_da3211_1_preprocessed.mat');

t = t(:);
N = length(t);

% The script assumes:
% z_m    = N x 12 measurements
% u_IMU  = N x 6 IMU inputs
% x_true = N x 18 true states
% Q      = 6 x 6 IMU-noise covariance
% R      = 12 x 12 measurement-noise covariance

%% Run the normal EKF

[x_est, P_est, z_est] = runEKF(t, z_m, u_IMU, Q, R, g);

%% Run the EKF again with airspeed noise sigma = 5 m/s

z_high = z_m;
V_true = sqrt(x_true(:,4).^2 + x_true(:,5).^2 + x_true(:,6).^2);
z_high(:,10) = V_true + 5*randn(N,1);

R_high = R;
R_high(10,10) = 5^2;

[x_est_high, P_est_high, z_est_high] = ...
    runEKF(t, z_high, u_IMU, Q, R_high, g);

%% Errors and RMSE

error_normal = x_est - x_true;
error_high   = x_est_high - x_true;

error_normal(:,7:9) = wrapAngle(error_normal(:,7:9));
error_high(:,7:9)   = wrapAngle(error_high(:,7:9));

RMSE_normal = sqrt(mean(error_normal.^2,1));
RMSE_high   = sqrt(mean(error_high.^2,1));

fprintf('\nRMSE of the first 9 states\n');
fprintf('State       Normal EKF       High airspeed noise\n');
for i = 1:9
    fprintf('%2d          %10.4f          %10.4f\n', ...
        i, RMSE_normal(i), RMSE_high(i));
end

%% Plot the first 9 states

stateName = {'x_E','y_E','z_E','u','v','w','phi','theta','psi'};

figure;
for i = 1:9
    subplot(3,3,i);
    hold on;
    grid on;

    if i >= 7
        plot(t, x_true(:,i)*180/pi, 'b');
        plot(t, x_est(:,i)*180/pi, 'r--');
        plot(t, x_est_high(:,i)*180/pi, 'k:');
        ylabel('deg');
    else
        plot(t, x_true(:,i), 'b');
        plot(t, x_est(:,i), 'r--');
        plot(t, x_est_high(:,i), 'k:');
    end

    title(stateName{i});
    xlabel('Time [s]');
end
legend('True','Normal EKF','High airspeed noise');

%% Plot estimated biases and wind

estimateName = {'lambda Ax','lambda Ay','lambda Az', ...
                'lambda p','lambda q','lambda r', ...
                'W_x','W_y','W_z'};

figure;
for i = 1:9
    stateNumber = i + 9;
    subplot(3,3,i);
    hold on;
    grid on;

    if stateNumber >= 13 && stateNumber <= 15
        plot(t, x_est(:,stateNumber)*180/pi, 'r');
        plot(t, x_est_high(:,stateNumber)*180/pi, 'k--');
        ylabel('deg/s');
    else
        plot(t, x_est(:,stateNumber), 'r');
        plot(t, x_est_high(:,stateNumber), 'k--');
    end

    title(estimateName{i});
    xlabel('Time [s]');
end
legend('Normal EKF','High airspeed noise');

%% Plot measured and estimated airspeed

figure;
hold on;
grid on;
plot(t, z_m(:,10), 'k.');
plot(t, z_est(:,10), 'r');
plot(t, z_high(:,10), 'b.');
plot(t, z_est_high(:,10), 'g');
xlabel('Time [s]');
ylabel('Airspeed [m/s]');
legend('Normal measurement','Normal EKF', ...
       'Noisy measurement','Noisy EKF');

%% Plot one-standard-deviation bounds for the normal EKF

sigma = zeros(N,18);
for k = 1:N
    sigma(k,:) = sqrt(diag(P_est(:,:,k))).';
end

figure;
for i = 1:9
    subplot(3,3,i);
    hold on;
    grid on;

    plot(t, error_normal(:,i), 'b');
    plot(t, sigma(:,i), 'r--');
    plot(t, -sigma(:,i), 'r--');
    title(stateName{i});
    xlabel('Time [s]');
end
legend('Error','+1 sigma','-1 sigma');


%% Extended Kalman Filter

function [x_history, P_history, z_history] = ...
    runEKF(t, z, u_IMU, Q, R, g)

    N  = length(t);
    nx = 18;

    % Initial state estimate
    x = zeros(nx,1);
    x(1:3) = z(1,1:3).';

    V     = z(1,10);
    alpha = z(1,11);
    beta  = z(1,12);

    x(4) = V*cos(alpha)*cos(beta);
    x(5) = V*sin(beta);
    x(6) = V*sin(alpha)*cos(beta);
    x(7:9) = z(1,7:9).';

    % Initial state uncertainty
    d2r = pi/180;
    initial_std = [10 10 10, 5 5 5, ...
                   2*d2r 2*d2r 2*d2r, ...
                   0.1 0.1 0.1, ...
                   0.05*d2r 0.05*d2r 0.05*d2r, ...
                   20 20 20];

    P = diag(initial_std.^2);

    x_history = zeros(N,nx);
    P_history = zeros(nx,nx,N);
    z_history = zeros(N,12);

    x_history(1,:) = x.';
    P_history(:,:,1) = P;
    z_history(1,:) = measurementModel(x).';

    for k = 1:N-1
        dt = t(k+1) - t(k);
        input = u_IMU(k,:).';

        % 1. Predict the state using Euler integration
        x_pred = x + dynamics(x,input,g)*dt;
        x_pred(7:9) = wrapAngle(x_pred(7:9));

        % 2. Calculate the Jacobians
        F = jacobian(@(state) dynamics(state,input,g), x_pred);
        H = jacobian(@measurementModel, x_pred);

        % 3. Discretize the linearized model
        G = noiseMatrix(x_pred);
        Phi = eye(nx) + F*dt;
        Gamma = G*dt;

        % 4. Predict the covariance
        P_pred = Phi*P*Phi.' + Gamma*Q*Gamma.';

        % Predict the measurement
        z_pred = measurementModel(x_pred);
        innovation = z(k+1,:).' - z_pred;
        innovation([7 8 9 11 12]) = ...
            wrapAngle(innovation([7 8 9 11 12]));

        S = H*P_pred*H.' + R;

        % 5. Calculate the Kalman gain
        K = (P_pred*H.')/S;

        % 6. Correct the state estimate
        x = x_pred + K*innovation;
        x(7:9) = wrapAngle(x(7:9));

        % 7. Correct the covariance
        P = (eye(nx) - K*H)*P_pred;
        P = (P + P.')/2;

        x_history(k+1,:) = x.';
        P_history(:,:,k+1) = P;
        z_history(k+1,:) = measurementModel(x).';
    end
end


%% Nonlinear aircraft equations

function xdot = dynamics(x,um,g)

    u = x(4);  v = x(5);  w = x(6);
    phi = x(7); theta = x(8); psi = x(9);

    % Remove estimated IMU biases
    Ax = um(1) - x(10);
    Ay = um(2) - x(11);
    Az = um(3) - x(12);
    p  = um(4) - x(13);
    q  = um(5) - x(14);
    r  = um(6) - x(15);

    Wx = x(16); Wy = x(17); Wz = x(18);

    cph = cos(phi);   sph = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cps = cos(psi);   sps = sin(psi);

    xdot = zeros(18,1);

    % Position rates
    xdot(1) = (u*cth + (v*sph+w*cph)*sth)*cps ...
              - (v*cph-w*sph)*sps + Wx;
    xdot(2) = (u*cth + (v*sph+w*cph)*sth)*sps ...
              + (v*cph-w*sph)*cps + Wy;
    xdot(3) = -u*sth + (v*sph+w*cph)*cth + Wz;

    % Body velocity rates
    xdot(4) = Ax - g*sth + r*v - q*w;
    xdot(5) = Ay + g*sph*cth + p*w - r*u;
    xdot(6) = Az + g*cph*cth + q*u - p*v;

    % Euler-angle rates
    xdot(7) = p + q*sph*tan(theta) + r*cph*tan(theta);
    xdot(8) = q*cph - r*sph;
    xdot(9) = (q*sph + r*cph)/cth;

    % Biases and wind are assumed constant
    xdot(10:18) = 0;
end


%% Measurement model

function z = measurementModel(x)

    u = x(4);  v = x(5);  w = x(6);
    phi = x(7); theta = x(8); psi = x(9);

    cph = cos(phi);   sph = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cps = cos(psi);   sps = sin(psi);

    z = zeros(12,1);

    % GPS position
    z(1:3) = x(1:3);

    % GPS velocity
    z(4) = (u*cth + (v*sph+w*cph)*sth)*cps ...
           - (v*cph-w*sph)*sps + x(16);
    z(5) = (u*cth + (v*sph+w*cph)*sth)*sps ...
           + (v*cph-w*sph)*cps + x(17);
    z(6) = -u*sth + (v*sph+w*cph)*cth + x(18);

    % Attitude
    z(7:9) = x(7:9);

    % Airspeed, angle of attack and sideslip angle
    z(10) = sqrt(u^2 + v^2 + w^2);
    z(11) = atan2(w,u);
    z(12) = atan2(v,sqrt(u^2+w^2));
end


%% Numerical Jacobian

function J = jacobian(fun,x)

    y = fun(x);
    J = zeros(length(y),length(x));
    step = 1e-6;

    for i = 1:length(x)
        x2 = x;
        x2(i) = x2(i) + step;
        J(:,i) = (fun(x2)-y)/step;
    end
end


%% IMU-noise input matrix

function G = noiseMatrix(x)

    u = x(4); v = x(5); w = x(6);
    phi = x(7); theta = x(8);

    G = zeros(18,6);

    G(4,1) = 1;
    G(5,2) = 1;
    G(6,3) = 1;

    G(5,4) =  w;
    G(6,4) = -v;
    G(7,4) =  1;

    G(4,5) = -w;
    G(6,5) =  u;
    G(7,5) =  sin(phi)*tan(theta);
    G(8,5) =  cos(phi);
    G(9,5) =  sin(phi)/cos(theta);

    G(4,6) =  v;
    G(5,6) = -u;
    G(7,6) =  cos(phi)*tan(theta);
    G(8,6) = -sin(phi);
    G(9,6) =  cos(phi)/cos(theta);
end


%% Keep angles between -pi and pi

function angle = wrapAngle(angle)
    angle = mod(angle + pi,2*pi) - pi;
end
