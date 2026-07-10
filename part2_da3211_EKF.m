clc;
clear;
close all;

%% PART 2 - EXTENDED KALMAN FILTER FOR MANEUVER DA3211_1
%
% This script follows the seven EKF steps shown in the assignment:
%
%   1. One-step-ahead state prediction
%   2. Calculation of the state and measurement Jacobians
%   3. Discretization of the continuous linearized model
%   4. Prediction of the state-error covariance matrix
%   5. Calculation of the Kalman gain
%   6. Measurement correction
%   7. Correction of the state-error covariance matrix
%
% The variable names deliberately follow the mathematical notation:
%
%   x_hat_k_k       = x-hat_{k,k}
%   x_hat_k1_k      = x-hat_{k+1,k}
%   x_hat_k1_k1     = x-hat_{k+1,k+1}
%
%   P_k_k           = P_{k,k}
%   P_k1_k          = P_{k+1,k}
%   P_k1_k1         = P_{k+1,k+1}
%
%   F_x_k           = F_x evaluated during the interval k -> k+1
%   H_x_k1          = H_x evaluated at x-hat_{k+1,k}
%   Phi_k1_k        = Phi_{k+1,k}
%   Gamma_k1_k      = Gamma_{k+1,k}
%   K_k1            = K_{k+1}
%
% State vector:
%
% x = [x_E, y_E, z_E, u, v, w, phi, theta, psi, ...
%      lambda_Ax, lambda_Ay, lambda_Az, ...
%      lambda_p, lambda_q, lambda_r, W_x, W_y, W_z]^T
%
% IMU input vector:
%
% u_m = [A_x,m, A_y,m, A_z,m, p_m, q_m, r_m]^T
%
% Measurement vector:
%
% z = [x_GPS, y_GPS, z_GPS, Vx_GPS, Vy_GPS, Vz_GPS, ...
%      phi_GPS, theta_GPS, psi_GPS, V, alpha, beta]^T

rng(2026);  % Keep the high-noise comparison reproducible.

%% 1. Load and check the Part 1 data

dataFile = 'part1_da3211_1_preprocessed.mat';
data = load(dataFile);

t           = data.t(:);
z_m         = data.z_m;
u_IMU       = data.u_IMU;
Q_d         = data.Q;
R           = data.R;
x_true      = data.x_true;
lambda_true = data.lambda_true(:);
wind_true   = data.wind_true(:);

N = length(t);

%% 2. Define the two airspeed-noise cases

nominalCase.name    = 'Nominal airspeed noise';
nominalCase.z       = z_m;
nominalCase.R       = R;
nominalCase.sigma_V = sqrt(R(10, 10));

highNoiseCase.name    = 'High airspeed noise: sigma_V = 5 m/s';
highNoiseCase.z       = z_m;
highNoiseCase.R       = R;
highNoiseCase.sigma_V = 5.0;

% Make only measurement 10, the true airspeed measurement, noisier.
V_true = sqrt(sum(x_true(:, 4:6).^2, 2));

highNoiseCase.z(:, 10) = ...
    V_true + highNoiseCase.sigma_V * randn(N, 1);

highNoiseCase.R(10, 10) = highNoiseCase.sigma_V^2;

%% 3. Run the EKF for both cases

fprintf('\n============================================================\n');
fprintf('Running the nominal EKF case\n');
fprintf('============================================================\n');

resultsNominal = runEKFCase( ...
    t, nominalCase.z, u_IMU, Q_d, nominalCase.R, x_true);

fprintf('\n============================================================\n');
fprintf('Running the high-airspeed-noise EKF case\n');
fprintf('============================================================\n');

resultsHighV = runEKFCase( ...
    t, highNoiseCase.z, u_IMU, Q_d, highNoiseCase.R, x_true);

%% 4. Print a compact comparison

stateNames = { ...
    'x_E'; 'y_E'; 'z_E'; ...
    'u'; 'v'; 'w'; ...
    'phi'; 'theta'; 'psi'; ...
    'lambda_Ax'; 'lambda_Ay'; 'lambda_Az'; ...
    'lambda_p'; 'lambda_q'; 'lambda_r'; ...
    'W_x'; 'W_y'; 'W_z'};

stateUnits = { ...
    'm'; 'm'; 'm'; ...
    'm/s'; 'm/s'; 'm/s'; ...
    'deg'; 'deg'; 'deg'; ...
    'm/s^2'; 'm/s^2'; 'm/s^2'; ...
    'deg/s'; 'deg/s'; 'deg/s'; ...
    'm/s'; 'm/s'; 'm/s'};

displayScale = ones(18, 1);
displayScale(7:9)   = 180 / pi;
displayScale(13:15) = 180 / pi;

trueFinal = x_true(end, :).'.* displayScale;

nominalFinal = resultsNominal.x_hat_cor(end, :).'.* displayScale;
highVFinal   = resultsHighV.x_hat_cor(end, :).'.* displayScale;

nominalFinalError = resultsNominal.final_error_cor .* displayScale;
highVFinalError   = resultsHighV.final_error_cor   .* displayScale;

nominalRMSE = resultsNominal.RMSE_cor .* displayScale;
highVRMSE   = resultsHighV.RMSE_cor   .* displayScale;

comparisonTable = table( ...
    stateNames, stateUnits, trueFinal, nominalFinal, highVFinal, ...
    nominalFinalError, highVFinalError, nominalRMSE, highVRMSE, ...
    'VariableNames', { ...
        'State', 'Unit', 'TrueFinal', 'NominalFinal', 'HighVFinal', ...
        'NominalFinalError', 'HighVFinalError', ...
        'NominalRMSE', 'HighVRMSE'});

fprintf('\n===================== EKF COMPARISON =====================\n');
disp(comparisonTable);

useForNIS = t > (t(1) + 1);

fprintf('\nMean normalized innovation squared after the first second:\n');
fprintf('  Nominal case: %.4f\n', ...
    mean(resultsNominal.NIS(useForNIS), 'omitnan'));
fprintf('  High-V case:  %.4f\n', ...
    mean(resultsHighV.NIS(useForNIS), 'omitnan'));
fprintf('  Reference mean for 12 measurements: 12\n');

printBiasAndWindResults( ...
    resultsNominal, resultsHighV, lambda_true, wind_true);

%% 5. Save everything needed for later analysis

outputFile = 'part2_da3211_EKF_results_humanized.mat';

save(outputFile, ...
    't', 'x_true', 'u_IMU', 'Q_d', ...
    'nominalCase', 'highNoiseCase', ...
    'resultsNominal', 'resultsHighV', ...
    'comparisonTable');

fprintf('\nSaved the EKF results to: %s\n', outputFile);

%% 6. Create the assignment plots

plotMeasurementComparison( ...
    t, nominalCase.z, resultsNominal.z_hat_cor, ...
    'Nominal case: raw and EKF-filtered measurements');

plotMeasurementComparison( ...
    t, highNoiseCase.z, resultsHighV.z_hat_cor, ...
    'High-V-noise case: raw and EKF-filtered measurements');

plotMainStateComparison( ...
    t, x_true, resultsNominal.x_hat_cor, resultsHighV.x_hat_cor);

plotMainStateErrors( ...
    t, resultsNominal.error_cor, resultsHighV.error_cor);

plotBiasAndWindEstimates( ...
    t, resultsNominal.x_hat_cor, resultsHighV.x_hat_cor, ...
    lambda_true, wind_true);

plotSelectedStandardDeviations(t, resultsNominal.sigma_x_cor);

plotErrorsWithSigmaBounds( ...
    t, resultsNominal.error_cor, resultsNominal.sigma_x_cor);

plotNormalizedInnovations( ...
    t, resultsNominal.normalized_innovation, ...
    'Nominal case: normalized innovations');

plotNormalizedInnovations( ...
    t, resultsHighV.normalized_innovation, ...
    'High-V-noise case: normalized innovations');

plotNISComparison(t, resultsNominal.NIS, resultsHighV.NIS);


%% ========================================================================
%  EKF IMPLEMENTATION
%  ========================================================================

function results = runEKFCase(t, z, u_IMU, Q_d, R, x_true)
%RUNEKFCASE Run one complete 18-state extended Kalman filter.
%
% The loop is written as k -> k+1 so that the MATLAB variables match the
% symbols in the EKF diagram as closely as possible.

    n_x = 18;
    n_z = 12;
    N   = numel(t);

    z_k   = z.';
    u_m_k = u_IMU.';

    %% Initial state estimate x-hat_{0,0}

    x_hat_0_0 = zeros(n_x, 1);

    % Position comes directly from the first GPS position measurement.
    x_hat_0_0(1:3) = z(1, 1:3).';

    % Convert the first V-alpha-beta airdata sample into body-axis velocity.
    V_0     = z(1, 10);
    alpha_0 = z(1, 11);
    beta_0  = z(1, 12);

    x_hat_0_0(4) = V_0 * cos(alpha_0) * cos(beta_0);
    x_hat_0_0(5) = V_0 * sin(beta_0);
    x_hat_0_0(6) = V_0 * sin(alpha_0) * cos(beta_0);

    % Initial attitude comes from the first attitude measurement.
    x_hat_0_0(7:9) = z(1, 7:9).';

    % No prior information is assumed for the six IMU biases and wind.
    x_hat_0_0(10:18) = 0;

    %% Initial covariance P_{0,0}

    deg2rad = pi / 180;

    initialStandardDeviation = [ ...
        10; 10; 10; ...                        % position [m]
        5; 5; 5; ...                           % body velocity [m/s]
        2*deg2rad; 2*deg2rad; 2*deg2rad; ...   % Euler angles [rad]
        0.1; 0.1; 0.1; ...                     % accelerometer biases [m/s^2]
        0.05*deg2rad; 0.05*deg2rad; ...
        0.05*deg2rad; ...                       % gyro biases [rad/s]
        20; 20; 20];                           % wind components [m/s]

    P_0_0 = diag(initialStandardDeviation.^2);

    %% Preallocate all histories

    x_hat_pred = zeros(N, n_x);
    x_hat_cor  = zeros(N, n_x);

    P_pred = zeros(n_x, n_x, N);
    P_cor  = zeros(n_x, n_x, N);

    sigma_x_pred = zeros(n_x, N);
    sigma_x_cor  = zeros(n_x, N);

    z_hat_pred = zeros(n_z, N);
    z_hat_cor  = zeros(n_z, N);

    innovation            = zeros(n_z, N);
    normalized_innovation = zeros(n_z, N);
    innovation_covariance = zeros(n_z, n_z, N);
    NIS                    = zeros(1, N);

    %% Store the initial estimate

    x_hat_k_k = x_hat_0_0;
    P_k_k     = P_0_0;

    x_hat_pred(1, :) = x_hat_k_k.';
    x_hat_cor(1, :)  = x_hat_k_k.';

    P_pred(:, :, 1) = P_k_k;
    P_cor(:, :, 1)  = P_k_k;

    sigma_x_pred(:, 1) = sqrt(max(diag(P_k_k), 0));
    sigma_x_cor(:, 1)  = sqrt(max(diag(P_k_k), 0));

    z_hat_pred(:, 1) = measurementModel12(x_hat_k_k);
    z_hat_cor(:, 1)  = z_hat_pred(:, 1);

    innovation(:, 1) = makeInnovation(z_k(:, 1), z_hat_pred(:, 1));

    H_x_0 = numericalJacobian(@measurementModel12, x_hat_k_k);
    S_0   = H_x_0 * P_k_k * H_x_0.' + R;
    S_0   = makeSymmetric(S_0);

    innovation_covariance(:, :, 1) = S_0;
    normalized_innovation(:, 1) = ...
        innovation(:, 1) ./ sqrt(max(diag(S_0), eps));
    NIS(1) = innovation(:, 1).' * (S_0 \ innovation(:, 1));

    %% Main EKF loop: propagate sample k to sample k+1

    tic;

    for k = 1:(N - 1)

        k1 = k + 1;
        dt = t(k1) - t(k);

        u_star_k = u_m_k(:, k);

        % ---------------------------------------------------------------
        % STEP 1: ONE-STEP-AHEAD STATE PREDICTION
        %
        % x-hat_{k+1,k} = x-hat_{k,k}
        %                 + integral f(x-hat, u*, t) dt
        % ---------------------------------------------------------------

        x_hat_k1_k = rk4Step( ...
            @navigationDynamics18, x_hat_k_k, u_star_k, t(k), dt);

        x_hat_k1_k(7:9) = wrapToPiLocal(x_hat_k1_k(7:9));

        % ---------------------------------------------------------------
        % STEP 2: CALCULATE THE JACOBIANS
        %
        % F_x = partial f / partial x
        % H_x = partial h / partial x
        % ---------------------------------------------------------------

        F_x_k = numericalJacobian( ...
            @(x) navigationDynamics18(t(k1), x, u_star_k), ...
            x_hat_k1_k);

        H_x_k1 = numericalJacobian( ...
            @measurementModel12, x_hat_k1_k);

        % G_k maps the six IMU noise components into the 18 state rates.
        G_k = imuNoiseMapping(x_hat_k1_k);

        % ---------------------------------------------------------------
        % STEP 3: DISCRETIZE THE LINEARIZED MODEL
        %
        % Phi_{k+1,k} and Gamma_{k+1,k}
        % ---------------------------------------------------------------

        [Phi_k1_k, Gamma_k1_k] = ...
            discretizeStateAndNoiseMatrices(F_x_k, G_k, dt);

        % ---------------------------------------------------------------
        % STEP 4: PREDICT THE STATE-ERROR COVARIANCE
        %
        % P_{k+1,k} = Phi P_{k,k} Phi^T
        %             + Gamma Q_{d,k} Gamma^T
        % ---------------------------------------------------------------

        P_k1_k = ...
            Phi_k1_k * P_k_k * Phi_k1_k.' + ...
            Gamma_k1_k * Q_d * Gamma_k1_k.';

        P_k1_k = makeSymmetric(P_k1_k);

        % Predict what the sensor vector should be at k+1.
        z_hat_k1_k = measurementModel12(x_hat_k1_k);

        % Innovation: actual measurement minus predicted measurement.
        innovation_k1 = makeInnovation(z_k(:, k1), z_hat_k1_k);

        % Innovation covariance:
        % S_{k+1} = H P_{k+1,k} H^T + R_{k+1}
        S_k1 = H_x_k1 * P_k1_k * H_x_k1.' + R;
        S_k1 = makeSymmetric(S_k1);

        % ---------------------------------------------------------------
        % STEP 5: CALCULATE THE KALMAN GAIN
        %
        % K_{k+1} = P_{k+1,k} H^T
        %           [H P_{k+1,k} H^T + R]^{-1}
        %
        % Matrix right division is used instead of forming inv(S_k1).
        % ---------------------------------------------------------------

        K_k1 = (P_k1_k * H_x_k1.') / S_k1;

        % ---------------------------------------------------------------
        % STEP 6: MEASUREMENT UPDATE
        %
        % x-hat_{k+1,k+1}
        % = x-hat_{k+1,k} + K_{k+1}[z_{k+1} - h(x-hat_{k+1,k})]
        % ---------------------------------------------------------------

        x_hat_k1_k1 = x_hat_k1_k + K_k1 * innovation_k1;
        x_hat_k1_k1(7:9) = wrapToPiLocal(x_hat_k1_k1(7:9));

        % ---------------------------------------------------------------
        % STEP 7: CORRECT THE STATE-ERROR COVARIANCE
        %
        % The diagram shows:
        % P_{k+1,k+1} = [I - K H]P_{k+1,k}.
        %
        % The Joseph form below is algebraically equivalent in exact
        % arithmetic, but is more robust to numerical round-off.
        % ---------------------------------------------------------------

        I_x = eye(n_x);

        P_k1_k1 = ...
            (I_x - K_k1 * H_x_k1) * P_k1_k * ...
            (I_x - K_k1 * H_x_k1).' + ...
            K_k1 * R * K_k1.';

        P_k1_k1 = makeSymmetric(P_k1_k1);

        %% Store everything before moving to the next sample

        x_hat_pred(k1, :) = x_hat_k1_k.';
        x_hat_cor(k1, :)  = x_hat_k1_k1.';

        P_pred(:, :, k1) = P_k1_k;
        P_cor(:, :, k1)  = P_k1_k1;

        sigma_x_pred(:, k1) = sqrt(max(diag(P_k1_k), 0));
        sigma_x_cor(:, k1)  = sqrt(max(diag(P_k1_k1), 0));

        z_hat_pred(:, k1) = z_hat_k1_k;
        z_hat_cor(:, k1)  = measurementModel12(x_hat_k1_k1);

        innovation(:, k1) = innovation_k1;
        innovation_covariance(:, :, k1) = S_k1;

        normalized_innovation(:, k1) = ...
            innovation_k1 ./ sqrt(max(diag(S_k1), eps));

        NIS(k1) = innovation_k1.' * (S_k1 \ innovation_k1);

        % The corrected estimate becomes x-hat_{k,k} for the next pass.
        x_hat_k_k = x_hat_k1_k1;
        P_k_k     = P_k1_k1;

        if mod(k1, 1000) == 0 || k1 == N
            fprintf('  k = %d of %d, time = %.2f s\n', k1, N, t(k1));
        end
    end

    elapsedTime = toc;
    fprintf('EKF completed in %.2f seconds.\n', elapsedTime);

    %% Estimation errors and RMSE

    error_pred = x_hat_pred - x_true;
    error_cor  = x_hat_cor  - x_true;

    error_pred(:, 7:9) = wrapToPiLocal(error_pred(:, 7:9));
    error_cor(:, 7:9)  = wrapToPiLocal(error_cor(:, 7:9));

    final_error_pred = x_hat_pred(end, :).' - x_true(end, :).';
    final_error_cor  = x_hat_cor(end, :).'  - x_true(end, :).';

    final_error_pred(7:9) = wrapToPiLocal(final_error_pred(7:9));
    final_error_cor(7:9)  = wrapToPiLocal(final_error_cor(7:9));

    RMSE_pred = sqrt(mean(error_pred.^2, 1)).';
    RMSE_cor  = sqrt(mean(error_cor.^2, 1)).';

    %% Return one structured result

    results.x_hat_pred = x_hat_pred;
    results.x_hat_cor  = x_hat_cor;

    results.P_pred = P_pred;
    results.P_cor  = P_cor;

    results.sigma_x_pred = sigma_x_pred;
    results.sigma_x_cor  = sigma_x_cor;

    results.z_hat_pred = z_hat_pred;
    results.z_hat_cor  = z_hat_cor;

    results.innovation            = innovation;
    results.normalized_innovation = normalized_innovation;
    results.innovation_covariance = innovation_covariance;
    results.NIS                    = NIS;

    results.error_pred = error_pred;
    results.error_cor  = error_cor;

    results.final_error_pred = final_error_pred;
    results.final_error_cor  = final_error_cor;

    results.RMSE_pred = RMSE_pred;
    results.RMSE_cor  = RMSE_cor;

    results.elapsedTime = elapsedTime;

    % Compatibility aliases for older versions of the Part 2 script.
    results.xhat_pred = x_hat_pred;
    results.xhat      = x_hat_cor;
    results.STD_x_pred = sigma_x_pred;
    results.STD_x_cor  = sigma_x_cor;
    results.Z_pred = z_hat_pred;
    results.Z_filt = z_hat_cor;
    results.Innov = innovation;
    results.NormInnov = normalized_innovation;
    results.EstErr_pred = error_pred;
    results.EstErr_cor  = error_cor;

end


%% ========================================================================
%  NONLINEAR PROCESS AND MEASUREMENT MODELS
%  ========================================================================

function x_dot = navigationDynamics18(~, x, u_m)
%NAVIGATIONDYNAMICS18 Continuous nonlinear 18-state navigation model.

    g = 9.80665;

    u     = x(4);
    v     = x(5);
    w     = x(6);

    phi   = x(7);
    theta = x(8);
    psi   = x(9);

    lambda_Ax = x(10);
    lambda_Ay = x(11);
    lambda_Az = x(12);

    lambda_p = x(13);
    lambda_q = x(14);
    lambda_r = x(15);

    W_x = x(16);
    W_y = x(17);
    W_z = x(18);

    A_x_m = u_m(1);
    A_y_m = u_m(2);
    A_z_m = u_m(3);

    p_m = u_m(4);
    q_m = u_m(5);
    r_m = u_m(6);

    % Remove the currently estimated IMU biases.
    A_x = A_x_m - lambda_Ax;
    A_y = A_y_m - lambda_Ay;
    A_z = A_z_m - lambda_Az;

    p = p_m - lambda_p;
    q = q_m - lambda_q;
    r = r_m - lambda_r;

    c_phi = cos(phi);
    s_phi = sin(phi);

    c_theta = cos(theta);
    s_theta = sin(theta);
    t_theta = tan(theta);

    c_psi = cos(psi);
    s_psi = sin(psi);

    x_dot = zeros(18, 1);

    % Position rates in the Earth/navigation frame.
    x_dot(1) = ...
        (u*c_theta + (v*s_phi + w*c_phi)*s_theta)*c_psi ...
        - (v*c_phi - w*s_phi)*s_psi + W_x;

    x_dot(2) = ...
        (u*c_theta + (v*s_phi + w*c_phi)*s_theta)*s_psi ...
        + (v*c_phi - w*s_phi)*c_psi + W_y;

    x_dot(3) = ...
        -u*s_theta + (v*s_phi + w*c_phi)*c_theta + W_z;

    % Body-axis velocity dynamics.
    x_dot(4) = A_x - g*s_theta       + r*v - q*w;
    x_dot(5) = A_y + g*s_phi*c_theta + p*w - r*u;
    x_dot(6) = A_z + g*c_phi*c_theta + q*u - p*v;

    % Euler-angle kinematics.
    x_dot(7) = p + q*s_phi*t_theta + r*c_phi*t_theta;
    x_dot(8) = q*c_phi - r*s_phi;
    x_dot(9) = (q*s_phi + r*c_phi) / c_theta;

    % Biases and wind are modeled as constant states.
    x_dot(10:18) = 0;

end


function z_hat = measurementModel12(x)
%MEASUREMENTMODEL12 Map the 18 states to the 12 sensor outputs.

    x_E = x(1);
    y_E = x(2);
    z_E = x(3);

    u = x(4);
    v = x(5);
    w = x(6);

    phi   = x(7);
    theta = x(8);
    psi   = x(9);

    W_x = x(16);
    W_y = x(17);
    W_z = x(18);

    c_phi = cos(phi);
    s_phi = sin(phi);

    c_theta = cos(theta);
    s_theta = sin(theta);

    c_psi = cos(psi);
    s_psi = sin(psi);

    z_hat = zeros(12, 1);

    % GPS position.
    z_hat(1:3) = [x_E; y_E; z_E];

    % GPS velocity in the Earth/navigation frame.
    z_hat(4) = ...
        (u*c_theta + (v*s_phi + w*c_phi)*s_theta)*c_psi ...
        - (v*c_phi - w*s_phi)*s_psi + W_x;

    z_hat(5) = ...
        (u*c_theta + (v*s_phi + w*c_phi)*s_theta)*s_psi ...
        + (v*c_phi - w*s_phi)*c_psi + W_y;

    z_hat(6) = ...
        -u*s_theta + (v*s_phi + w*c_phi)*c_theta + W_z;

    % Attitude measurements.
    z_hat(7:9) = [phi; theta; psi];

    % Airdata measurements.
    z_hat(10) = sqrt(u^2 + v^2 + w^2);
    z_hat(11) = atan2(w, u);
    z_hat(12) = atan2(v, sqrt(u^2 + w^2));

end


%% ========================================================================
%  EKF SUPPORT FUNCTIONS
%  ========================================================================

function x_next = rk4Step(f, x_now, u, t_now, dt)
%RK4STEP One classical fourth-order Runge-Kutta integration step.

    k_1 = f(t_now,          x_now,                u);
    k_2 = f(t_now + dt/2,   x_now + dt*k_1/2,     u);
    k_3 = f(t_now + dt/2,   x_now + dt*k_2/2,     u);
    k_4 = f(t_now + dt,     x_now + dt*k_3,       u);

    x_next = x_now + dt*(k_1 + 2*k_2 + 2*k_3 + k_4)/6;

end


function J = numericalJacobian(fun, x)
%NUMERICALJACOBIAN Central-difference Jacobian of a vector-valued function.

    y = fun(x);

    n_y = numel(y);
    n_x = numel(x);

    J = zeros(n_y, n_x);

    for iState = 1:n_x
        step = 1e-6 * max(1, abs(x(iState)));

        x_plus  = x;
        x_minus = x;

        x_plus(iState)  = x_plus(iState)  + step;
        x_minus(iState) = x_minus(iState) - step;

        y_plus  = fun(x_plus);
        y_minus = fun(x_minus);

        J(:, iState) = (y_plus - y_minus) / (2*step);
    end

end


function G = imuNoiseMapping(x)
%IMUNOISEMAPPING Map accelerometer and gyro noise into the state rates.

    u = x(4);
    v = x(5);
    w = x(6);

    phi   = x(7);
    theta = x(8);

    c_phi = cos(phi);
    s_phi = sin(phi);

    c_theta = cos(theta);
    t_theta = tan(theta);

    G = zeros(18, 6);

    % Accelerometer-noise channels.
    G(4, 1) = 1;
    G(5, 2) = 1;
    G(6, 3) = 1;

    % p-channel gyro noise.
    G(5, 4) =  w;
    G(6, 4) = -v;
    G(7, 4) =  1;

    % q-channel gyro noise.
    G(4, 5) = -w;
    G(6, 5) =  u;
    G(7, 5) =  s_phi*t_theta;
    G(8, 5) =  c_phi;
    G(9, 5) =  s_phi/c_theta;

    % r-channel gyro noise.
    G(4, 6) =  v;
    G(5, 6) = -u;
    G(7, 6) =  c_phi*t_theta;
    G(8, 6) = -s_phi;
    G(9, 6) =  c_phi/c_theta;

end


function [Phi, Gamma] = discretizeStateAndNoiseMatrices(F_x, G, dt)
%DISCRETIZESTATEANDNOISEMATRICES Exact ZOH discretization via expm.

    n_x = size(F_x, 1);
    n_w = size(G, 2);

    augmentedMatrix = [ ...
        F_x,          G; ...
        zeros(n_w, n_x), zeros(n_w, n_w)];

    discreteAugmentedMatrix = expm(augmentedMatrix * dt);

    Phi   = discreteAugmentedMatrix(1:n_x, 1:n_x);
    Gamma = discreteAugmentedMatrix(1:n_x, n_x + 1:n_x + n_w);

end


function nu = makeInnovation(z_measured, z_predicted)
%MAKEINNOVATION Compute z - z-hat and wrap angular residuals.

    nu = z_measured - z_predicted;

    % phi, theta, psi, alpha and beta are angular measurements.
    angularMeasurementIndices = [7, 8, 9, 11, 12];

    nu(angularMeasurementIndices) = ...
        wrapToPiLocal(nu(angularMeasurementIndices));

end


function angle = wrapToPiLocal(angle)
%WRAPTOPILOCAL Wrap angles to the interval [-pi, pi).

    angle = mod(angle + pi, 2*pi) - pi;

end


function A = makeSymmetric(A)
%MAKESYMMETRIC Remove tiny numerical asymmetry from a covariance matrix.

    A = 0.5 * (A + A.');

end


%% ========================================================================
%  DATA AND REPORTING HELPERS
%  ========================================================================

function dataFile = findPart1DataFile()
%FINDPART1DATAFILE Locate the Part 1 MAT-file in the current folder.

    preferredName = 'part1_da3211_1_preprocessed.mat';

    if isfile(preferredName)
        dataFile = preferredName;
        return;
    end

    candidates = dir('part1_da3211_1_preprocessed*.mat');

    if isempty(candidates)
        error([ ...
            'No Part 1 data file was found. Put the preprocessed MAT-file ' ...
            'in the same folder as this script.']);
    end

    [~, newestIndex] = max([candidates.datenum]);
    dataFile = candidates(newestIndex).name;

end


function X = orientSamplesByRow(X, N, variableName)
%ORIENTSAMPLESBYROW Make every row correspond to one time sample.

    if size(X, 1) == N
        return;
    end

    if size(X, 2) == N
        X = X.';
        return;
    end

    error('%s does not have one dimension equal to length(t).', variableName);

end


function printBiasAndWindResults( ...
    resultsNominal, resultsHighV, lambda_true, wind_true)
%PRINTBIASANDWINDRESULTS Print the final estimated nuisance states.

    rad2deg = 180 / pi;

    fprintf('\n================ FINAL ESTIMATED BIASES ================\n');

    fprintf('\nAccelerometer biases [m/s^2]\n');
    fprintf('                    lambda_Ax      lambda_Ay      lambda_Az\n');
    fprintf('  Nominal:       %12.6f %12.6f %12.6f\n', ...
        resultsNominal.x_hat_cor(end, 10:12));
    fprintf('  High-V noise:  %12.6f %12.6f %12.6f\n', ...
        resultsHighV.x_hat_cor(end, 10:12));

    fprintf('\nGyro biases [deg/s]\n');
    fprintf('                    lambda_p       lambda_q       lambda_r\n');
    fprintf('  Nominal:       %12.6f %12.6f %12.6f\n', ...
        resultsNominal.x_hat_cor(end, 13:15) * rad2deg);
    fprintf('  High-V noise:  %12.6f %12.6f %12.6f\n', ...
        resultsHighV.x_hat_cor(end, 13:15) * rad2deg);

    fprintf('\n================== FINAL ESTIMATED WIND =================\n');
    fprintf('                    W_x            W_y            W_z\n');
    fprintf('  Nominal:       %12.6f %12.6f %12.6f\n', ...
        resultsNominal.x_hat_cor(end, 16:18));
    fprintf('  High-V noise:  %12.6f %12.6f %12.6f\n', ...
        resultsHighV.x_hat_cor(end, 16:18));

    if numel(lambda_true) == 6
        fprintf('\nTrue bias vector:\n');
        disp(lambda_true.');

        fprintf('Nominal final bias error:\n');
        disp(resultsNominal.x_hat_cor(end, 10:15).' - lambda_true);

        fprintf('High-V final bias error:\n');
        disp(resultsHighV.x_hat_cor(end, 10:15).' - lambda_true);
    end

    if numel(wind_true) == 3
        fprintf('\nTrue wind vector:\n');
        disp(wind_true.');

        fprintf('Nominal final wind error:\n');
        disp(resultsNominal.x_hat_cor(end, 16:18).' - wind_true);

        fprintf('High-V final wind error:\n');
        disp(resultsHighV.x_hat_cor(end, 16:18).' - wind_true);
    end

end


%% ========================================================================
%  PLOTTING HELPERS
%  ========================================================================

function plotMeasurementComparison(t, z_measured, z_filtered, figureTitle)

    measurementLabels = { ...
        'x_E [m]', 'y_E [m]', 'z_E [m]', ...
        'V_{x,GPS} [m/s]', 'V_{y,GPS} [m/s]', 'V_{z,GPS} [m/s]', ...
        '\phi [deg]', '\theta [deg]', '\psi [deg]', ...
        'V [m/s]', '\alpha [deg]', '\beta [deg]'};

    angularMeasurements = [7, 8, 9, 11, 12];

    figure('Name', figureTitle);
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(4, 3, 'TileSpacing', 'compact');

    for iMeasurement = 1:12
        nexttile;
        hold on;
        grid on;

        measuredSignal = z_measured(:, iMeasurement);
        filteredSignal = z_filtered(iMeasurement, :).';

        if ismember(iMeasurement, angularMeasurements)
            measuredSignal = measuredSignal * 180/pi;
            filteredSignal = filteredSignal * 180/pi;
        end

        plot(t, measuredSignal, 'k.', 'MarkerSize', 4);
        plot(t, filteredSignal, 'r', 'LineWidth', 1.2);

        xlabel('Time [s]');
        ylabel(measurementLabels{iMeasurement});
        legend('Raw measurement', 'h(x-hat_{k,k})', ...
               'Location', 'best');
    end

    sgtitle(figureTitle);

end


function plotMainStateComparison(t, x_true, x_nominal, x_highV)

    stateLabels = { ...
        'x_E [m]', 'y_E [m]', 'z_E [m]', ...
        'u [m/s]', 'v [m/s]', 'w [m/s]', ...
        '\phi [deg]', '\theta [deg]', '\psi [deg]'};

    figure('Name', 'True state and corrected EKF estimates');
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for iState = 1:9
        nexttile;
        hold on;
        grid on;

        trueSignal    = x_true(:, iState);
        nominalSignal = x_nominal(:, iState);
        highVSignal   = x_highV(:, iState);

        if iState >= 7
            trueSignal    = trueSignal * 180/pi;
            nominalSignal = nominalSignal * 180/pi;
            highVSignal   = highVSignal * 180/pi;
        end

        plot(t, trueSignal, 'b', 'LineWidth', 1.3);
        plot(t, nominalSignal, 'r--', 'LineWidth', 1.2);
        plot(t, highVSignal, 'k:', 'LineWidth', 1.2);

        xlabel('Time [s]');
        ylabel(stateLabels{iState});
        legend('True', 'Nominal EKF', 'High-V-noise EKF', ...
               'Location', 'best');
    end

    sgtitle('True state versus corrected EKF estimates');

end


function plotMainStateErrors(t, errorNominal, errorHighV)

    stateLabels = { ...
        'x_E [m]', 'y_E [m]', 'z_E [m]', ...
        'u [m/s]', 'v [m/s]', 'w [m/s]', ...
        '\phi [deg]', '\theta [deg]', '\psi [deg]'};

    figure('Name', 'Corrected EKF state-estimation errors');
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for iState = 1:9
        nexttile;
        hold on;
        grid on;

        nominalSignal = errorNominal(:, iState);
        highVSignal   = errorHighV(:, iState);

        if iState >= 7
            nominalSignal = nominalSignal * 180/pi;
            highVSignal   = highVSignal * 180/pi;
        end

        plot(t, nominalSignal, 'r', 'LineWidth', 1.1);
        plot(t, highVSignal, 'k--', 'LineWidth', 1.1);

        xlabel('Time [s]');
        ylabel(['Error in ', stateLabels{iState}]);
        legend('Nominal', 'High-V-noise', 'Location', 'best');
    end

    sgtitle('Corrected EKF state-estimation errors');

end


function plotBiasAndWindEstimates( ...
    t, x_nominal, x_highV, lambda_true, wind_true)

    estimateLabels = { ...
        '\lambda_{A_x} [m/s^2]', ...
        '\lambda_{A_y} [m/s^2]', ...
        '\lambda_{A_z} [m/s^2]', ...
        '\lambda_p [deg/s]', ...
        '\lambda_q [deg/s]', ...
        '\lambda_r [deg/s]', ...
        'W_x [m/s]', 'W_y [m/s]', 'W_z [m/s]'};

    estimateIndices = 10:18;

    figure('Name', 'Estimated IMU biases and wind');
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for iEstimate = 1:9
        nexttile;
        hold on;
        grid on;

        stateIndex = estimateIndices(iEstimate);

        nominalSignal = x_nominal(:, stateIndex);
        highVSignal   = x_highV(:, stateIndex);

        isGyroBias = iEstimate >= 4 && iEstimate <= 6;

        if isGyroBias
            nominalSignal = nominalSignal * 180/pi;
            highVSignal   = highVSignal * 180/pi;
        end

        plot(t, nominalSignal, 'r', 'LineWidth', 1.2);
        plot(t, highVSignal, 'k--', 'LineWidth', 1.2);

        trueLineWasAdded = false;

        if numel(lambda_true) == 6 && iEstimate <= 6
            trueValue = lambda_true(iEstimate);

            if isGyroBias
                trueValue = trueValue * 180/pi;
            end

            yline(trueValue, 'b:', 'LineWidth', 1.3);
            trueLineWasAdded = true;

        elseif numel(wind_true) == 3 && iEstimate >= 7
            trueValue = wind_true(iEstimate - 6);
            yline(trueValue, 'b:', 'LineWidth', 1.3);
            trueLineWasAdded = true;
        end

        xlabel('Time [s]');
        ylabel(estimateLabels{iEstimate});

        if trueLineWasAdded
            legend('Nominal', 'High-V-noise', 'True', 'Location', 'best');
        else
            legend('Nominal', 'High-V-noise', 'Location', 'best');
        end
    end

    sgtitle('Estimated IMU biases and wind states');

end


function plotSelectedStandardDeviations(t, sigma_x)

    selectedStates = [1, 4, 7, 10, 11, 12, 16, 17, 18];

    selectedLabels = { ...
        '\sigma_{x_E} [m]', ...
        '\sigma_u [m/s]', ...
        '\sigma_\phi [deg]', ...
        '\sigma_{\lambda A_x} [m/s^2]', ...
        '\sigma_{\lambda A_y} [m/s^2]', ...
        '\sigma_{\lambda A_z} [m/s^2]', ...
        '\sigma_{W_x} [m/s]', ...
        '\sigma_{W_y} [m/s]', ...
        '\sigma_{W_z} [m/s]'};

    figure('Name', 'Nominal corrected state standard deviations');
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for iPlot = 1:numel(selectedStates)
        nexttile;
        hold on;
        grid on;

        stateIndex = selectedStates(iPlot);
        sigmaSignal = sigma_x(stateIndex, :).';

        if stateIndex == 7
            sigmaSignal = sigmaSignal * 180/pi;
        end

        plot(t, sigmaSignal, 'r', 'LineWidth', 1.2);

        xlabel('Time [s]');
        ylabel(selectedLabels{iPlot});
    end

    sgtitle('Nominal EKF corrected state standard deviations');

end


function plotErrorsWithSigmaBounds(t, error_cor, sigma_x_cor)

    selectedStates = [1, 4, 7, 10, 11, 12, 16, 17, 18];

    selectedLabels = { ...
        'x_E error [m]', ...
        'u error [m/s]', ...
        '\phi error [deg]', ...
        '\lambda_{A_x} error [m/s^2]', ...
        '\lambda_{A_y} error [m/s^2]', ...
        '\lambda_{A_z} error [m/s^2]', ...
        'W_x error [m/s]', ...
        'W_y error [m/s]', ...
        'W_z error [m/s]'};

    figure('Name', 'Nominal errors with one-sigma bounds');
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for iPlot = 1:numel(selectedStates)
        nexttile;
        hold on;
        grid on;

        stateIndex = selectedStates(iPlot);

        errorSignal = error_cor(:, stateIndex);
        sigmaSignal = sigma_x_cor(stateIndex, :).';

        if stateIndex == 7
            errorSignal = errorSignal * 180/pi;
            sigmaSignal = sigmaSignal * 180/pi;
        end

        plot(t, errorSignal, 'b', 'LineWidth', 1.0);
        plot(t, sigmaSignal, 'r--', 'LineWidth', 1.0);
        plot(t, -sigmaSignal, 'r--', 'LineWidth', 1.0);

        xlabel('Time [s]');
        ylabel(selectedLabels{iPlot});
        legend('Estimation error', '+1\sigma', '-1\sigma', ...
               'Location', 'best');
    end

    sgtitle('Nominal EKF errors and one-sigma covariance bounds');

end


function plotNormalizedInnovations(t, normalizedInnovation, figureTitle)

    figure('Name', figureTitle);
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(4, 3, 'TileSpacing', 'compact');

    for iMeasurement = 1:12
        nexttile;
        hold on;
        grid on;

        plot(t, normalizedInnovation(iMeasurement, :), ...
             'b', 'LineWidth', 1.0);

        yline(0, 'k-');
        yline(3, 'r:');
        yline(-3, 'r:');

        xlabel('Time [s]');
        ylabel(sprintf('\\nu_%d / \\sigma_{\\nu_%d}', ...
                       iMeasurement, iMeasurement));
    end

    sgtitle(figureTitle);

end


function plotNISComparison(t, NIS_nominal, NIS_highV)

    figure('Name', 'Normalized innovation squared comparison');
    set(gcf, 'Position', [100, 100, 1000, 500]);
    hold on;
    grid on;

    plot(t, NIS_nominal, 'r', 'LineWidth', 1.1);
    plot(t, NIS_highV, 'k--', 'LineWidth', 1.1);

    yline(12, 'b--', 'Expected mean: n_z = 12');
    yline(21.03, 'm:', 'Approximate 95% chi-square bound');
    yline(26.22, 'm:', 'Approximate 99% chi-square bound');

    xlabel('Time [s]');
    ylabel('NIS');
    legend('Nominal', 'High-V-noise', 'Expected mean', ...
           'Location', 'best');
    title('Normalized innovation squared comparison');

end
