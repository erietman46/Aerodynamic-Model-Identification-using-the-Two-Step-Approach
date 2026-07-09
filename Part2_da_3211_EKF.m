clc;
clear;
close all;

%% ========================================================================
%  PART 2 - EKF FOR DA3211_1 MANEUVER
%
%  This script:
%  1. Loads the Part 1 preprocessed file for da3211_1.
%  2. Runs the nominal EKF with sigma_V = 0.2 m/s.
%  3. Runs a second EKF case with sigma_V = 5 m/s.
%  4. Compares true state, EKF prediction, EKF correction, raw measurements,
%     filtered measurements, bias estimates, wind estimates, innovations,
%     NIS, and RMSE.
%
%  STATES:
%  x = [xE yE zE u v w phi theta psi ...
%       lambda_Ax lambda_Ay lambda_Az lambda_p lambda_q lambda_r ...
%       Wx Wy Wz]^T
%
%  INPUTS:
%  u_m = [Ax_m Ay_m Az_m p_m q_m r_m]^T
%
%  MEASUREMENTS:
%  z_m = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS ...
%         phi_GPS theta_GPS psi_GPS V_m alpha_m beta_m]^T
%% ========================================================================

rng(2026);

%% ========================= LOAD DATA ====================================

dataFile = 'part1_da3211_1_preprocessed.mat';

if exist(dataFile, 'file') ~= 2
    error('Could not find part1_da3211_1_preprocessed.mat in the current folder.');
end

data = load(dataFile);

t       = data.t(:);
z_m     = data.z_m;
u_IMU   = data.u_IMU;
Q       = data.Q;
R       = data.R;
x_true  = data.x_true;

if isfield(data, 'lambda_true')
    lambda_true = data.lambda_true(:);
else
    lambda_true = [];
end

if isfield(data, 'wind_true')
    wind_true = data.wind_true(:);
else
    wind_true = [];
end

N = length(t);

% Ensure correct dimensions
if size(z_m,1) ~= N && size(z_m,2) == N
    z_m = z_m.';
end

if size(u_IMU,1) ~= N && size(u_IMU,2) == N
    u_IMU = u_IMU.';
end

if size(x_true,1) ~= N && size(x_true,2) == N
    x_true = x_true.';
end

%% ========================= CASE DEFINITIONS ==============================

% Nominal case
caseNominal.name = 'Nominal airspeed noise';
caseNominal.z_m  = z_m;
caseNominal.R    = R;
caseNominal.sigma_V = sqrt(R(10,10));

% High airspeed-noise case required by the assignment
caseHighV.name = 'High airspeed noise, sigma_V = 5 m/s';
caseHighV.z_m  = z_m;
caseHighV.R    = R;
caseHighV.sigma_V = 5.0;

% Replace only the true airspeed measurement with a noisier one
V_true = sqrt(x_true(:,4).^2 + x_true(:,5).^2 + x_true(:,6).^2);

caseHighV.z_m(:,10) = V_true + caseHighV.sigma_V*randn(N,1);
caseHighV.R(10,10)  = caseHighV.sigma_V^2;

%% ========================= RUN BOTH EKF CASES ============================

fprintf('\n============================================================\n');
fprintf('Running nominal EKF case...\n');
fprintf('============================================================\n');

resultsNominal = runEKFcase(t, caseNominal.z_m, u_IMU, Q, caseNominal.R, x_true);

fprintf('\n============================================================\n');
fprintf('Running high-V-noise EKF case...\n');
fprintf('============================================================\n');

resultsHighV = runEKFcase(t, caseHighV.z_m, u_IMU, Q, caseHighV.R, x_true);

%% ========================= PRINT RESULTS =================================

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

scale = ones(18,1);
scale(7:9)   = 180/pi;
scale(13:15) = 180/pi;

RMSE_nominal = resultsNominal.RMSE_cor .* scale;
RMSE_highV   = resultsHighV.RMSE_cor   .* scale;

Final_nominal = resultsNominal.xhat(end,:).' .* scale;
Final_highV   = resultsHighV.xhat(end,:).'   .* scale;
True_final    = x_true(end,:).'              .* scale;

FinalErr_nominal = resultsNominal.final_cor_error .* scale;
FinalErr_highV   = resultsHighV.final_cor_error   .* scale;

comparisonTable = table( ...
    stateNames, ...
    stateUnits, ...
    True_final, ...
    Final_nominal, ...
    Final_highV, ...
    FinalErr_nominal, ...
    FinalErr_highV, ...
    RMSE_nominal, ...
    RMSE_highV, ...
    'VariableNames', { ...
    'State', ...
    'Unit', ...
    'TrueFinal', ...
    'EKFNominalFinal', ...
    'EKFHighVFinal', ...
    'NominalFinalError', ...
    'HighVFinalError', ...
    'NominalRMSE', ...
    'HighVRMSE'});

fprintf('\n================ EKF COMPARISON TABLE ================\n');
disp(comparisonTable);

fprintf('\nMean NIS after first second:\n');
fprintf('Nominal case: %.4f\n', mean(resultsNominal.NIS(t > 1), 'omitnan'));
fprintf('High-V case:  %.4f\n', mean(resultsHighV.NIS(t > 1), 'omitnan'));
fprintf('Expected mean NIS is approximately n_z = 12.\n');

fprintf('\n================ FINAL ESTIMATED BIASES ================\n');

fprintf('\nNominal accelerometer biases [m/s^2]:\n');
fprintf('lambda_Ax = %.8f\n', resultsNominal.xhat(end,10));
fprintf('lambda_Ay = %.8f\n', resultsNominal.xhat(end,11));
fprintf('lambda_Az = %.8f\n', resultsNominal.xhat(end,12));

fprintf('\nHigh-V accelerometer biases [m/s^2]:\n');
fprintf('lambda_Ax = %.8f\n', resultsHighV.xhat(end,10));
fprintf('lambda_Ay = %.8f\n', resultsHighV.xhat(end,11));
fprintf('lambda_Az = %.8f\n', resultsHighV.xhat(end,12));

fprintf('\nNominal gyro biases [deg/s]:\n');
fprintf('lambda_p = %.8f\n', resultsNominal.xhat(end,13)*180/pi);
fprintf('lambda_q = %.8f\n', resultsNominal.xhat(end,14)*180/pi);
fprintf('lambda_r = %.8f\n', resultsNominal.xhat(end,15)*180/pi);

fprintf('\nHigh-V gyro biases [deg/s]:\n');
fprintf('lambda_p = %.8f\n', resultsHighV.xhat(end,13)*180/pi);
fprintf('lambda_q = %.8f\n', resultsHighV.xhat(end,14)*180/pi);
fprintf('lambda_r = %.8f\n', resultsHighV.xhat(end,15)*180/pi);

fprintf('\n================ FINAL ESTIMATED WIND ================\n');

fprintf('\nNominal wind [m/s]:\n');
fprintf('Wx = %.8f\n', resultsNominal.xhat(end,16));
fprintf('Wy = %.8f\n', resultsNominal.xhat(end,17));
fprintf('Wz = %.8f\n', resultsNominal.xhat(end,18));

fprintf('\nHigh-V wind [m/s]:\n');
fprintf('Wx = %.8f\n', resultsHighV.xhat(end,16));
fprintf('Wy = %.8f\n', resultsHighV.xhat(end,17));
fprintf('Wz = %.8f\n', resultsHighV.xhat(end,18));

if ~isempty(lambda_true)
    fprintf('\nTrue IMU biases:\n');
    disp(lambda_true.');

    fprintf('\nNominal final bias error:\n');
    disp(resultsNominal.xhat(end,10:15).' - lambda_true);

    fprintf('\nHigh-V final bias error:\n');
    disp(resultsHighV.xhat(end,10:15).' - lambda_true);
end

if ~isempty(wind_true)
    fprintf('\nTrue wind:\n');
    disp(wind_true.');

    fprintf('\nNominal final wind error:\n');
    disp(resultsNominal.xhat(end,16:18).' - wind_true);

    fprintf('\nHigh-V final wind error:\n');
    disp(resultsHighV.xhat(end,16:18).' - wind_true);
end

%% ========================= SAVE RESULTS ==================================

save('part2_da3211_EKF_results_corrected.mat', ...
    't', ...
    'x_true', ...
    'u_IMU', ...
    'caseNominal', ...
    'caseHighV', ...
    'resultsNominal', ...
    'resultsHighV', ...
    'comparisonTable');

fprintf('\nSaved corrected EKF results to part2_da3211_EKF_results_corrected.mat\n');

%% ========================= PLOTTING ======================================

measurementNames = { ...
    'x_E [m]', ...
    'y_E [m]', ...
    'z_E [m]', ...
    'u_{GPS} [m/s]', ...
    'v_{GPS} [m/s]', ...
    'w_{GPS} [m/s]', ...
    '\phi [deg]', ...
    '\theta [deg]', ...
    '\psi [deg]', ...
    'V [m/s]', ...
    '\alpha [deg]', ...
    '\beta [deg]'};

angleMeasurementIndex = [7 8 9 11 12];

stateNames9 = { ...
    'x_E [m]', ...
    'y_E [m]', ...
    'z_E [m]', ...
    'u [m/s]', ...
    'v [m/s]', ...
    'w [m/s]', ...
    '\phi [deg]', ...
    '\theta [deg]', ...
    '\psi [deg]'};

%% 1. Raw measurements versus filtered measurements, nominal case

figure('Name','Nominal raw measurements vs filtered measurements');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(4,3,'TileSpacing','compact');

for i = 1:12
    nexttile;
    hold on;
    grid on;

    measuredSignal = caseNominal.z_m(:,i);
    filteredSignal = resultsNominal.Z_filt(i,:).';

    if ismember(i, angleMeasurementIndex)
        measuredSignal = measuredSignal*180/pi;
        filteredSignal = filteredSignal*180/pi;
    end

    plot(t, measuredSignal, 'k.', 'MarkerSize', 4);
    plot(t, filteredSignal, 'r', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(measurementNames{i});
    legend('Raw measurement', 'Filtered measurement h(x_{k|k})', 'Location', 'best');
end

sgtitle('Nominal case: raw measurements versus filtered EKF measurements');

%% 2. Raw measurements versus filtered measurements, high V-noise case

figure('Name','High-V raw measurements vs filtered measurements');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(4,3,'TileSpacing','compact');

for i = 1:12
    nexttile;
    hold on;
    grid on;

    measuredSignal = caseHighV.z_m(:,i);
    filteredSignal = resultsHighV.Z_filt(i,:).';

    if ismember(i, angleMeasurementIndex)
        measuredSignal = measuredSignal*180/pi;
        filteredSignal = filteredSignal*180/pi;
    end

    plot(t, measuredSignal, 'k.', 'MarkerSize', 4);
    plot(t, filteredSignal, 'r', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(measurementNames{i});
    legend('Raw measurement', 'Filtered measurement h(x_{k|k})', 'Location', 'best');
end

sgtitle('High-V-noise case: raw measurements versus filtered EKF measurements');

%% 3. True state versus nominal and high-V corrected estimates

figure('Name','True state vs nominal/high-V EKF estimates');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:9
    nexttile;
    hold on;
    grid on;

    trueSignal = x_true(:,i);
    nominalSignal = resultsNominal.xhat(:,i);
    highVSignal = resultsHighV.xhat(:,i);

    if i >= 7
        trueSignal = trueSignal*180/pi;
        nominalSignal = nominalSignal*180/pi;
        highVSignal = highVSignal*180/pi;
    end

    plot(t, trueSignal, 'b', 'LineWidth', 1.3);
    plot(t, nominalSignal, 'r--', 'LineWidth', 1.2);
    plot(t, highVSignal, 'k:', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(stateNames9{i});
    legend('True', 'Nominal EKF', 'High-V-noise EKF', 'Location', 'best');
end

sgtitle('True state versus EKF corrected estimates');

%% 4. Estimation error comparison for main states

figure('Name','State estimation errors: nominal vs high-V');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:9
    nexttile;
    hold on;
    grid on;

    errNominal = resultsNominal.EstErr_cor(:,i);
    errHighV   = resultsHighV.EstErr_cor(:,i);

    if i >= 7
        errNominal = errNominal*180/pi;
        errHighV   = errHighV*180/pi;
    end

    plot(t, errNominal, 'r', 'LineWidth', 1.1);
    plot(t, errHighV, 'k--', 'LineWidth', 1.1);

    xlabel('Time [s]');
    ylabel(['Error in ', stateNames9{i}]);
    legend('Nominal', 'High-V-noise', 'Location', 'best');
end

sgtitle('Corrected EKF state estimation errors');

%% 5. Bias and wind estimates, nominal versus high-V

estimateNames = { ...
    '\lambda_{A_x} [m/s^2]', ...
    '\lambda_{A_y} [m/s^2]', ...
    '\lambda_{A_z} [m/s^2]', ...
    '\lambda_p [deg/s]', ...
    '\lambda_q [deg/s]', ...
    '\lambda_r [deg/s]', ...
    'W_x [m/s]', ...
    'W_y [m/s]', ...
    'W_z [m/s]'};

estimateIndex = 10:18;

figure('Name','Estimated biases and wind: nominal vs high-V');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(3,3,'TileSpacing','compact');

for j = 1:9
    nexttile;
    hold on;
    grid on;

    idx = estimateIndex(j);

    nominalSignal = resultsNominal.xhat(:,idx);
    highVSignal   = resultsHighV.xhat(:,idx);

    if j >= 4 && j <= 6
        nominalSignal = nominalSignal*180/pi;
        highVSignal   = highVSignal*180/pi;
    end

    plot(t, nominalSignal, 'r', 'LineWidth', 1.2);
    plot(t, highVSignal, 'k--', 'LineWidth', 1.2);

    if ~isempty(lambda_true) && j <= 6
        trueValue = lambda_true(j);
        if j >= 4
            trueValue = trueValue*180/pi;
        end
        yline(trueValue, 'b:', 'LineWidth', 1.3);
        legend('Nominal', 'High-V-noise', 'True', 'Location', 'best');

    elseif ~isempty(wind_true) && j >= 7
        trueValue = wind_true(j-6);
        yline(trueValue, 'b:', 'LineWidth', 1.3);
        legend('Nominal', 'High-V-noise', 'True', 'Location', 'best');

    else
        legend('Nominal', 'High-V-noise', 'Location', 'best');
    end

    xlabel('Time [s]');
    ylabel(estimateNames{j});
end

sgtitle('Estimated IMU biases and wind states');

%% 6. State standard deviations, nominal case

selectedStates = [1 4 7 10 11 12 16 17 18];

selectedNames = { ...
    '\sigma_{x_E} [m]', ...
    '\sigma_u [m/s]', ...
    '\sigma_\phi [deg]', ...
    '\sigma_{\lambda A_x} [m/s^2]', ...
    '\sigma_{\lambda A_y} [m/s^2]', ...
    '\sigma_{\lambda A_z} [m/s^2]', ...
    '\sigma_{W_x} [m/s]', ...
    '\sigma_{W_y} [m/s]', ...
    '\sigma_{W_z} [m/s]'};

figure('Name','Nominal state standard deviations');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(3,3,'TileSpacing','compact');

for j = 1:length(selectedStates)
    nexttile;
    hold on;
    grid on;

    idx = selectedStates(j);
    stdSignal = resultsNominal.STD_x_cor(idx,:);

    if idx == 7
        stdSignal = stdSignal*180/pi;
    end

    plot(t, stdSignal, 'r', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(selectedNames{j});
end

sgtitle('Nominal EKF corrected state standard deviations');

%% 7. Estimation errors with +/- 1 sigma, nominal case

figure('Name','Nominal estimation errors with 1 sigma bounds');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(3,3,'TileSpacing','compact');

for j = 1:length(selectedStates)
    nexttile;
    hold on;
    grid on;

    idx = selectedStates(j);

    errSignal = resultsNominal.EstErr_cor(:,idx);
    stdSignal = resultsNominal.STD_x_cor(idx,:).';

    if idx == 7
        errSignal = errSignal*180/pi;
        stdSignal = stdSignal*180/pi;
    end

    plot(t, errSignal, 'b', 'LineWidth', 1.0);
    plot(t, stdSignal, 'r--', 'LineWidth', 1.0);
    plot(t, -stdSignal, 'r--', 'LineWidth', 1.0);

    xlabel('Time [s]');
    ylabel(selectedNames{j});
    legend('Error', '+1\sigma', '-1\sigma', 'Location', 'best');
end

sgtitle('Nominal EKF corrected estimation error with covariance bounds');

%% 8. Proper normalized innovations, nominal case

figure('Name','Nominal normalized innovations');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(4,3,'TileSpacing','compact');

for i = 1:12
    nexttile;
    hold on;
    grid on;

    plot(t, resultsNominal.NormInnov(i,:), 'b', 'LineWidth', 1.0);
    yline(0, 'k-');
    yline(3, 'r:');
    yline(-3, 'r:');

    xlabel('Time [s]');
    ylabel(['\nu_', num2str(i), ' / \sigma_{\nu}']);
end

sgtitle('Nominal case: normalized innovation sequence');

%% 9. Proper normalized innovations, high-V case

figure('Name','High-V normalized innovations');
set(gcf, 'Position', [100 100 1300 850]);
tiledlayout(4,3,'TileSpacing','compact');

for i = 1:12
    nexttile;
    hold on;
    grid on;

    plot(t, resultsHighV.NormInnov(i,:), 'b', 'LineWidth', 1.0);
    yline(0, 'k-');
    yline(3, 'r:');
    yline(-3, 'r:');

    xlabel('Time [s]');
    ylabel(['\nu_', num2str(i), ' / \sigma_{\nu}']);
end

sgtitle('High-V-noise case: normalized innovation sequence');

%% 10. NIS comparison

figure('Name','NIS comparison');
set(gcf, 'Position', [100 100 1000 500]);
hold on;
grid on;

plot(t, resultsNominal.NIS, 'r', 'LineWidth', 1.1);
plot(t, resultsHighV.NIS, 'k--', 'LineWidth', 1.1);

yline(12, 'b--', 'Expected mean, n_z = 12');
yline(21.03, 'm:', 'Approx. 95% chi-square bound');
yline(26.22, 'm:', 'Approx. 99% chi-square bound');

xlabel('Time [s]');
ylabel('NIS');
legend('Nominal', 'High-V-noise', 'Expected mean', 'Location', 'best');
title('Normalized Innovation Squared comparison');

%% ========================================================================
%  LOCAL EKF FUNCTION
%% ========================================================================

function results = runEKFcase(t, z_m, u_IMU, Q_imu, R_meas, x_true)

    nx = 18;
    nz = 12;
    N  = length(t);
    d2r = pi/180;

    Z_k = z_m.';
    U_k = u_IMU.';

    %% ---------------------- Initialization ------------------------------

    x0 = zeros(nx,1);

    % Initial position from first GPS measurement
    x0(1) = z_m(1,1);
    x0(2) = z_m(1,2);
    x0(3) = z_m(1,3);

    % Initial body velocity from first airdata measurement
    V0     = z_m(1,10);
    alpha0 = z_m(1,11);
    beta0  = z_m(1,12);

    x0(4) = V0*cos(alpha0)*cos(beta0);
    x0(5) = V0*sin(beta0);
    x0(6) = V0*sin(alpha0)*cos(beta0);

    % Initial attitude from first GPS attitude measurement
    x0(7) = z_m(1,7);
    x0(8) = z_m(1,8);
    x0(9) = z_m(1,9);

    % Initial bias guesses
    x0(10:15) = 0;

    % Initial wind guesses
    x0(16:18) = 0;

    % Initial covariance
    std0 = [ ...
        10; 10; 10; ...                      % position uncertainty [m]
        5; 5; 5; ...                         % body velocity uncertainty [m/s]
        2*d2r; 2*d2r; 2*d2r; ...             % attitude uncertainty [rad]
        0.1; 0.1; 0.1; ...                   % accel bias uncertainty [m/s^2]
        0.05*d2r; 0.05*d2r; 0.05*d2r; ...    % gyro bias uncertainty [rad/s]
        20; 20; 20];                         % wind uncertainty [m/s]

    P0 = diag(std0.^2);

    %% ---------------------- Storage --------------------------------------

    xhat_pred = zeros(N,nx);
    xhat      = zeros(N,nx);

    P_pred = zeros(nx,nx,N);
    P_cor  = zeros(nx,nx,N);

    STD_x_pred = zeros(nx,N);
    STD_x_cor  = zeros(nx,N);

    Z_pred = zeros(nz,N);
    Z_filt = zeros(nz,N);

    Innov     = zeros(nz,N);
    NormInnov = zeros(nz,N);
    NIS       = zeros(1,N);

    %% ---------------------- Initial sample -------------------------------

    x_k_k = x0;
    P_k_k = P0;

    xhat_pred(1,:) = x_k_k.';
    xhat(1,:)      = x_k_k.';

    P_pred(:,:,1) = P_k_k;
    P_cor(:,:,1)  = P_k_k;

    STD_x_pred(:,1) = sqrt(diag(P_k_k));
    STD_x_cor(:,1)  = sqrt(diag(P_k_k));

    Z_pred(:,1) = measurementModel12_local(x_k_k);
    Z_filt(:,1) = measurementModel12_local(x_k_k);

    Innov(:,1) = innovation12_local(Z_k(:,1), Z_pred(:,1));

    H0 = numericalJacobian_local(@measurementModel12_local, x_k_k);
    S0 = H0*P_k_k*H0.' + R_meas;
    NormInnov(:,1) = Innov(:,1) ./ sqrt(diag(S0));
    NIS(1) = Innov(:,1).' * (S0 \ Innov(:,1));

    %% ---------------------- EKF loop -------------------------------------

    tic;

    for k = 2:N

        dt = t(k) - t(k-1);

        %% Prediction

        [~, x_k1_k] = rk4_local(@navDynamics18_local, x_k_k, U_k(:,k-1), [t(k-1), t(k)]);

        x_k1_k(7:9) = wrapPi_local(x_k1_k(7:9));

        Fx = numericalJacobian_local(@(x) navDynamics18_local(0, x, U_k(:,k-1)), x_k1_k);
        G  = imuNoiseMapping_local(x_k1_k);

        [Phi, Gamma] = c2d_local(Fx, G, dt);

        P_k1_k = Phi*P_k_k*Phi.' + Gamma*Q_imu*Gamma.';
        P_k1_k = 0.5*(P_k1_k + P_k1_k.');

        %% Correction

        Hx = numericalJacobian_local(@measurementModel12_local, x_k1_k);

        z_k1_k = measurementModel12_local(x_k1_k);
        innov  = innovation12_local(Z_k(:,k), z_k1_k);

        S = Hx*P_k1_k*Hx.' + R_meas;
        S = 0.5*(S + S.');

        K = P_k1_k*Hx.' / S;

        x_k1_k1 = x_k1_k + K*innov;
        x_k1_k1(7:9) = wrapPi_local(x_k1_k1(7:9));

        I = eye(nx);

        % Joseph-form covariance update
        P_k1_k1 = (I - K*Hx)*P_k1_k*(I - K*Hx).' + K*R_meas*K.';
        P_k1_k1 = 0.5*(P_k1_k1 + P_k1_k1.');

        %% Store

        xhat_pred(k,:) = x_k1_k.';
        xhat(k,:)      = x_k1_k1.';

        P_pred(:,:,k) = P_k1_k;
        P_cor(:,:,k)  = P_k1_k1;

        STD_x_pred(:,k) = sqrt(diag(P_k1_k));
        STD_x_cor(:,k)  = sqrt(diag(P_k1_k1));

        Z_pred(:,k) = z_k1_k;
        Z_filt(:,k) = measurementModel12_local(x_k1_k1);

        Innov(:,k) = innov;

        NormInnov(:,k) = innov ./ sqrt(diag(S));

        NIS(k) = innov.' * (S \ innov);

        x_k_k = x_k1_k1;
        P_k_k = P_k1_k1;

        if mod(k,1000) == 0 || k == N
            fprintf('EKF running: k = %d / %d, t = %.2f s\n', k, N, t(k));
        end
    end

    elapsedTime = toc;
    fprintf('EKF completed in %.2f seconds.\n', elapsedTime);

    %% ---------------------- Error analysis -------------------------------

    EstErr_pred = xhat_pred - x_true;
    EstErr_cor  = xhat      - x_true;

    EstErr_pred(:,7:9) = wrapPi_local(EstErr_pred(:,7:9));
    EstErr_cor(:,7:9)  = wrapPi_local(EstErr_cor(:,7:9));

    final_pred_error = xhat_pred(end,:).' - x_true(end,:).';
    final_cor_error  = xhat(end,:).'      - x_true(end,:).';

    final_pred_error(7:9) = wrapPi_local(final_pred_error(7:9));
    final_cor_error(7:9)  = wrapPi_local(final_cor_error(7:9));

    RMSE_pred = sqrt(mean(EstErr_pred.^2,1)).';
    RMSE_cor  = sqrt(mean(EstErr_cor.^2,1)).';

    %% ---------------------- Output structure -----------------------------

    results.xhat_pred = xhat_pred;
    results.xhat      = xhat;

    results.P_pred = P_pred;
    results.P_cor  = P_cor;

    results.STD_x_pred = STD_x_pred;
    results.STD_x_cor  = STD_x_cor;

    results.Z_pred = Z_pred;
    results.Z_filt = Z_filt;

    results.Innov     = Innov;
    results.NormInnov = NormInnov;
    results.NIS       = NIS;

    results.EstErr_pred = EstErr_pred;
    results.EstErr_cor  = EstErr_cor;

    results.final_pred_error = final_pred_error;
    results.final_cor_error  = final_cor_error;

    results.RMSE_pred = RMSE_pred;
    results.RMSE_cor  = RMSE_cor;

    results.elapsedTime = elapsedTime;

end

%% ========================================================================
%  LOCAL DYNAMICS AND MEASUREMENT FUNCTIONS
%% ========================================================================

function [tout, xout] = rk4_local(fhandle, x0, u, tspan)

    t0 = tspan(1);
    t1 = tspan(2);
    dt = t1 - t0;

    k1 = fhandle(t0,        x0,           u);
    k2 = fhandle(t0+dt/2,   x0+dt*k1/2,   u);
    k3 = fhandle(t0+dt/2,   x0+dt*k2/2,   u);
    k4 = fhandle(t1,        x0+dt*k3,     u);

    xout = x0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    tout = t1;

end

function xdot = navDynamics18_local(~, x, u_m)

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

    Wx = x(16);
    Wy = x(17);
    Wz = x(18);

    Ax_m = u_m(1);
    Ay_m = u_m(2);
    Az_m = u_m(3);

    p_m = u_m(4);
    q_m = u_m(5);
    r_m = u_m(6);

    Ax = Ax_m - lambda_Ax;
    Ay = Ay_m - lambda_Ay;
    Az = Az_m - lambda_Az;

    p = p_m - lambda_p;
    q = q_m - lambda_q;
    r = r_m - lambda_r;

    cphi = cos(phi);
    sphi = sin(phi);

    ctheta = cos(theta);
    stheta = sin(theta);
    ttheta = tan(theta);

    cpsi = cos(psi);
    spsi = sin(psi);

    xdot = zeros(18,1);

    % Position dynamics
    xdot(1) = (u*ctheta + (v*sphi + w*cphi)*stheta)*cpsi ...
            - (v*cphi - w*sphi)*spsi ...
            + Wx;

    xdot(2) = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
            + (v*cphi - w*sphi)*cpsi ...
            + Wy;

    xdot(3) = -u*stheta ...
            + (v*sphi + w*cphi)*ctheta ...
            + Wz;

    % Body velocity dynamics
    xdot(4) = Ax - g*stheta + r*v - q*w;
    xdot(5) = Ay + g*sphi*ctheta + p*w - r*u;
    xdot(6) = Az + g*cphi*ctheta + q*u - p*v;

    % Euler angle dynamics
    xdot(7) = p + q*sphi*ttheta + r*cphi*ttheta;
    xdot(8) = q*cphi - r*sphi;
    xdot(9) = (q*sphi + r*cphi)/ctheta;

    % Biases and wind are modeled as constants
    xdot(10:18) = 0;

end

function zhat = measurementModel12_local(x)

    xE = x(1);
    yE = x(2);
    zE = x(3);

    u = x(4);
    v = x(5);
    w = x(6);

    phi = x(7);
    theta = x(8);
    psi = x(9);

    Wx = x(16);
    Wy = x(17);
    Wz = x(18);

    cphi = cos(phi);
    sphi = sin(phi);

    ctheta = cos(theta);
    stheta = sin(theta);

    cpsi = cos(psi);
    spsi = sin(psi);

    zhat = zeros(12,1);

    % GPS position
    zhat(1) = xE;
    zhat(2) = yE;
    zhat(3) = zE;

    % GPS ground speed in navigation frame
    zhat(4) = (u*ctheta + (v*sphi + w*cphi)*stheta)*cpsi ...
            - (v*cphi - w*sphi)*spsi ...
            + Wx;

    zhat(5) = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
            + (v*cphi - w*sphi)*cpsi ...
            + Wy;

    zhat(6) = -u*stheta ...
            + (v*sphi + w*cphi)*ctheta ...
            + Wz;

    % GPS attitude
    zhat(7) = phi;
    zhat(8) = theta;
    zhat(9) = psi;

    % Airdata
    zhat(10) = sqrt(u^2 + v^2 + w^2);
    zhat(11) = atan2(w, u);
    zhat(12) = atan2(v, sqrt(u^2 + w^2));

end

function J = numericalJacobian_local(fun, x)

    y0 = fun(x);

    ny = length(y0);
    nx = length(x);

    J = zeros(ny,nx);

    for i = 1:nx
        dx = 1e-6*max(1,abs(x(i)));

        xp = x;
        xp(i) = xp(i) + dx;

        yp = fun(xp);

        J(:,i) = (yp - y0)/dx;
    end

end

function G = imuNoiseMapping_local(x)

    u = x(4);
    v = x(5);
    w = x(6);

    phi = x(7);
    theta = x(8);

    cphi = cos(phi);
    sphi = sin(phi);

    ctheta = cos(theta);
    ttheta = tan(theta);

    G = zeros(18,6);

    % Accelerometer noise mapping
    G(4,1) = 1;
    G(5,2) = 1;
    G(6,3) = 1;

    % Gyro p noise mapping
    G(5,4) = w;
    G(6,4) = -v;
    G(7,4) = 1;

    % Gyro q noise mapping
    G(4,5) = -w;
    G(6,5) = u;
    G(7,5) = sphi*ttheta;
    G(8,5) = cphi;
    G(9,5) = sphi/ctheta;

    % Gyro r noise mapping
    G(4,6) = v;
    G(5,6) = -u;
    G(7,6) = cphi*ttheta;
    G(8,6) = -sphi;
    G(9,6) = cphi/ctheta;

end

function [Phi, Gamma] = c2d_local(F, G, dt)

    n = size(F,1);
    m = size(G,2);

    A = [F, G;
         zeros(m,n), zeros(m,m)];

    B = expm(A*dt);

    Phi   = B(1:n,1:n);
    Gamma = B(1:n,n+1:n+m);

end

function nu = innovation12_local(z, zhat)

    nu = z - zhat;

    % Wrap angular residuals:
    % phi, theta, psi, alpha, beta
    angleIndex = [7 8 9 11 12];

    nu(angleIndex) = wrapPi_local(nu(angleIndex));

end

function angle = wrapPi_local(angle)

    angle = mod(angle + pi, 2*pi) - pi;

end