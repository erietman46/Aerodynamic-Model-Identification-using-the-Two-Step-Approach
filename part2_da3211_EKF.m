clc;
clear;
close all;

%% PART 2 - IMPROVED EXTENDED KALMAN FILTER FOR MANEUVER DA3211_1
%
% This submission-ready script implements and evaluates an 18-state EKF.
% It addresses all Part 2 requirements:
%
%   2.1  Complete EKF prediction and correction procedure.
%   2.2  Validation of accelerometer biases, gyro biases and wind states.
%   2.3  Comparison of raw sensor measurements with EKF-reconstructed data.
%   2.4  Paired comparison with true-airspeed noise increased to 5 m/s.
%   2.5  Convergence, innovation, NIS, NEES and whiteness diagnostics.
%
% State vector:
% x = [x_E, y_E, z_E, u, v, w, phi, theta, psi, ...
%      lambda_Ax, lambda_Ay, lambda_Az, ...
%      lambda_p, lambda_q, lambda_r, W_x, W_y, W_z]^T
%
% IMU input vector:
% u_m = [A_x,m, A_y,m, A_z,m, p_m, q_m, r_m]^T
%
% Measurement vector:
% z = [x_GPS, y_GPS, z_GPS, V_N,GPS, V_E,GPS, V_D,GPS, ...
%      phi_GPS, theta_GPS, psi_GPS, V, alpha, beta]^T

rng(2026, 'twister');

%% 0. User settings

settings.g = 9.8030;
% The supplied synthetic data are dynamically most consistent with
% g approximately equal to 9.803 m/s^2. The gravity diagnostic printed
% below independently checks this value using the known simulated truth.
% Do not silently tune g: report the value used and its effect on lambda_Az.

settings.burnInSeconds          = 1.0;
settings.confidenceLevel        = 0.95;
settings.maxInnovationLag       = 25;
settings.saveFigures            = true;
settings.figureResolution       = 200;
settings.usePairedAirspeedNoise = true;
settings.highAirspeedSigma      = 5.0;

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

outputFolder = fullfile(scriptFolder, 'part2_output');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% 1. Load and validate the Part 1 data

dataFile = findPart1DataFile(scriptFolder);
fprintf('Loading Part 1 data from:\n  %s\n', dataFile);

data = load(dataFile);
[data, t, z_m, u_IMU, Q_imu, R, x_true, lambda_true, wind_true] = ...
    validateAndExtractPart1Data(data);

N = numel(t);
fprintf('Number of samples: %d\n', N);
fprintf('Mean sample time: %.6f s\n', mean(diff(t)));

stateNames = { ...
    'x_E'; 'y_E'; 'z_E'; 'u'; 'v'; 'w'; ...
    'phi'; 'theta'; 'psi'; ...
    'lambda_Ax'; 'lambda_Ay'; 'lambda_Az'; ...
    'lambda_p'; 'lambda_q'; 'lambda_r'; ...
    'W_x'; 'W_y'; 'W_z'};

stateUnits = { ...
    'm'; 'm'; 'm'; 'm/s'; 'm/s'; 'm/s'; ...
    'deg'; 'deg'; 'deg'; ...
    'm/s^2'; 'm/s^2'; 'm/s^2'; ...
    'deg/s'; 'deg/s'; 'deg/s'; ...
    'm/s'; 'm/s'; 'm/s'};

measurementNames = { ...
    'x_GPS'; 'y_GPS'; 'z_GPS'; ...
    'V_N_GPS'; 'V_E_GPS'; 'V_D_GPS'; ...
    'phi_GPS'; 'theta_GPS'; 'psi_GPS'; ...
    'V_airdata'; 'alpha_airdata'; 'beta_airdata'};

%% 2. Gravity-model diagnostic

if numel(lambda_true) == 6
    gravityCheck = estimateEffectiveGravity( ...
        t, x_true, u_IMU, lambda_true);

    fprintf('\n================ GRAVITY MODEL CHECK ================\n');
    fprintf('Gravity used by EKF:                     %.6f m/s^2\n', settings.g);
    fprintf('Median gravity inferred from z-dynamics: %.6f m/s^2\n', ...
        gravityCheck.gFromZMedian);
    fprintf('Robust spread of inferred gravity:       %.6f m/s^2\n', ...
        gravityCheck.gFromZRobustStd);

    if abs(settings.g - gravityCheck.gFromZMedian) > 1.0e-3
        warning([ ...
            'The selected gravity differs from the value inferred from the ' ...
            'synthetic trajectory. This difference is expected to appear ' ...
            'mainly in the estimated vertical accelerometer bias.']);
    end
else
    gravityCheck = struct();
    warning('lambda_true is unavailable; the gravity diagnostic was skipped.');
end

%% 3. Define nominal and high-airspeed-noise cases

nominalCase.name    = 'Nominal airspeed noise';
nominalCase.z       = z_m;
nominalCase.R       = R;
nominalCase.sigma_V = sqrt(R(10, 10));

highNoiseCase.name    = sprintf( ...
    'High airspeed noise: sigma_V = %.1f m/s', settings.highAirspeedSigma);
highNoiseCase.z       = z_m;
highNoiseCase.R       = R;
highNoiseCase.sigma_V = settings.highAirspeedSigma;

V_true = sqrt(sum(x_true(:, 4:6).^2, 2));

if settings.usePairedAirspeedNoise && nominalCase.sigma_V > 0
    % Use the same standardized noise realization in both cases. This makes
    % the comparison fair: only the airspeed noise magnitude changes.
    standardizedVNoise = ...
        (nominalCase.z(:, 10) - V_true) / nominalCase.sigma_V;
    highNoiseCase.z(:, 10) = ...
        V_true + highNoiseCase.sigma_V * standardizedVNoise;
else
    highNoiseCase.z(:, 10) = ...
        V_true + highNoiseCase.sigma_V * randn(N, 1);
end

highNoiseCase.R(10, 10) = highNoiseCase.sigma_V^2;

%% 4. Run the EKF for both cases

fprintf('\n============================================================\n');
fprintf('Running nominal EKF case\n');
fprintf('============================================================\n');

resultsNominal = runEKFCase( ...
    t, nominalCase.z, u_IMU, Q_imu, nominalCase.R, x_true, settings);

fprintf('\n============================================================\n');
fprintf('Running high-airspeed-noise EKF case\n');
fprintf('============================================================\n');

resultsHighV = runEKFCase( ...
    t, highNoiseCase.z, u_IMU, Q_imu, highNoiseCase.R, x_true, settings);

%% 5. Quantitative assignment results

comparisonTable = buildStateComparisonTable( ...
    x_true, resultsNominal, resultsHighV, stateNames, stateUnits);

stateConsistencyNominal = buildStateConsistencyTable( ...
    t, resultsNominal, stateNames, stateUnits, settings);
stateConsistencyHighV = buildStateConsistencyTable( ...
    t, resultsHighV, stateNames, stateUnits, settings);

nuisanceTable = buildNuisanceConsistencyTable( ...
    resultsNominal, resultsHighV, lambda_true, wind_true);

innovationTableNominal = buildInnovationStatistics( ...
    t, resultsNominal.normalized_innovation, measurementNames, settings);
innovationTableHighV = buildInnovationStatistics( ...
    t, resultsHighV.normalized_innovation, measurementNames, settings);

consistencyNominal = analyseFilterConsistency( ...
    t, resultsNominal, settings);
consistencyHighV = analyseFilterConsistency( ...
    t, resultsHighV, settings);

fprintf('\n===================== STATE COMPARISON =====================\n');
disp(comparisonTable);

fprintf('\n========== BIAS AND WIND FINAL CONSISTENCY CHECK ==========\n');
disp(nuisanceTable);

fprintf('\n========== NOMINAL NORMALIZED-INNOVATION STATISTICS =======\n');
disp(innovationTableNominal);

fprintf('\n========== HIGH-V NORMALIZED-INNOVATION STATISTICS ========\n');
disp(innovationTableHighV);

printConsistencySummary('Nominal', consistencyNominal);
printConsistencySummary('High-V-noise', consistencyHighV);

fprintf('\n===== NOMINAL STATE-WISE COVARIANCE DIAGNOSTICS =====\n');
fprintf(['MeanSquaredNormalizedError should be of order 1 when each state ' ...
    'variance is well calibrated. Large values identify the states that ' ...
    'dominate the full-state NEES.\n']);
stateConsistencyNominalSorted = sortrows( ...
    stateConsistencyNominal, 'MeanSquaredNormalizedError', 'descend');
disp(stateConsistencyNominalSorted);

fprintf('\n===== HIGH-V STATE-WISE COVARIANCE DIAGNOSTICS =====\n');
stateConsistencyHighVSorted = sortrows( ...
    stateConsistencyHighV, 'MeanSquaredNormalizedError', 'descend');
disp(stateConsistencyHighVSorted);

%% 6. Create all assignment plots

plotMeasurementComparison( ...
    t, nominalCase.z, resultsNominal.z_hat_cor, ...
    'Nominal case: raw and EKF-reconstructed measurements');

plotMeasurementComparison( ...
    t, highNoiseCase.z, resultsHighV.z_hat_cor, ...
    'High-V-noise case: raw and EKF-reconstructed measurements');

plotMainStateComparison( ...
    t, x_true, resultsNominal.x_hat_cor, resultsHighV.x_hat_cor);

plotMainStateErrors( ...
    t, resultsNominal.error_cor, resultsHighV.error_cor);

plotBiasAndWindEstimates( ...
    t, resultsNominal.x_hat_cor, resultsHighV.x_hat_cor, ...
    lambda_true, wind_true);

plotNuisanceErrorsWith3Sigma( ...
    t, resultsNominal, lambda_true, wind_true, ...
    'Nominal bias and wind errors with 3-sigma bounds');

plotNuisanceErrorsWith3Sigma( ...
    t, resultsHighV, lambda_true, wind_true, ...
    'High-V-noise bias and wind errors with 3-sigma bounds');

plotSelectedStandardDeviations(t, resultsNominal.sigma_x_cor);

plotErrorsWithSigmaBounds( ...
    t, resultsNominal.error_cor, resultsNominal.sigma_x_cor);

plotNormalizedInnovations( ...
    t, resultsNominal.normalized_innovation, ...
    'Nominal case: normalized innovations');

plotNormalizedInnovations( ...
    t, resultsHighV.normalized_innovation, ...
    'High-V-noise case: normalized innovations');

plotInnovationAutocorrelations( ...
    t, resultsNominal.normalized_innovation, measurementNames, ...
    settings, ...
    'Nominal normalized-innovation autocorrelation');

plotNISComparison( ...
    t, resultsNominal.NIS, resultsHighV.NIS, ...
    consistencyNominal, consistencyHighV);

plotNEESComparison( ...
    t, resultsNominal.NEES, resultsHighV.NEES, ...
    consistencyNominal, consistencyHighV);

plotRMSEChange(comparisonTable);

%% 7. Save results, tables and figures

outputFile = fullfile(outputFolder, 'part2_da3211_EKF_improved_results.mat');

save(outputFile, ...
    'settings', 'gravityCheck', 't', 'x_true', 'u_IMU', 'Q_imu', ...
    'nominalCase', 'highNoiseCase', ...
    'resultsNominal', 'resultsHighV', ...
    'comparisonTable', 'stateConsistencyNominal', ...
    'stateConsistencyHighV', 'nuisanceTable', ...
    'innovationTableNominal', 'innovationTableHighV', ...
    'consistencyNominal', 'consistencyHighV');

writetable(comparisonTable, ...
    fullfile(outputFolder, 'state_comparison.csv'));
writetable(nuisanceTable, ...
    fullfile(outputFolder, 'bias_wind_consistency.csv'));
writetable(innovationTableNominal, ...
    fullfile(outputFolder, 'innovation_statistics_nominal.csv'));
writetable(innovationTableHighV, ...
    fullfile(outputFolder, 'innovation_statistics_highV.csv'));
writetable(stateConsistencyNominal, ...
    fullfile(outputFolder, 'state_consistency_nominal.csv'));
writetable(stateConsistencyHighV, ...
    fullfile(outputFolder, 'state_consistency_highV.csv'));

if settings.saveFigures
    saveAllFigures(outputFolder, settings.figureResolution);
end

fprintf('\nSaved improved Part 2 results to:\n  %s\n', outputFile);
fprintf('Tables and figures are in:\n  %s\n', outputFolder);

%% ========================================================================
%  EKF IMPLEMENTATION
%  ========================================================================

function results = runEKFCase(t, z, u_IMU, Q_d, R, x_true, settings)
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
            @(time, state, input) navigationDynamics18( ...
                time, state, input, settings.g), ...
            x_hat_k_k, u_star_k, t(k), dt);

        x_hat_k1_k(7:9) = wrapToPiLocal(x_hat_k1_k(7:9));

        % ---------------------------------------------------------------
        % STEP 2: CALCULATE THE JACOBIANS
        %
        % F_x = partial f / partial x
        % H_x = partial h / partial x
        % ---------------------------------------------------------------

        F_x_k = numericalJacobian( ...
            @(x) navigationDynamics18(t(k1), x, u_star_k, settings.g), ...
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

    % Normalized state errors and NEES are possible because this is a
    % simulated assignment for which x_true is known.
    normalized_state_error = ...
        error_cor ./ max(sigma_x_cor.', sqrt(eps));

    NEES = nan(1, N);
    for k = 1:N
        NEES(k) = safeQuadraticForm( ...
            error_cor(k, :).', P_cor(:, :, k));
    end

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

    results.normalized_state_error = normalized_state_error;
    results.NEES = NEES;

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

function x_dot = navigationDynamics18(~, x, u_m, g)
%NAVIGATIONDYNAMICS18 Continuous nonlinear 18-state navigation model.
% g is passed explicitly so the gravity value is documented and testable.

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

function dataFile = findPart1DataFile(scriptFolder)
%FINDPART1DATAFILE Locate the preprocessed Part 1 MAT-file portably.

    searchFolders = {scriptFolder, pwd};
    preferredName = 'part1_da3211_1_preprocessed.mat';

    for iFolder = 1:numel(searchFolders)
        candidate = fullfile(searchFolders{iFolder}, preferredName);
        if isfile(candidate)
            dataFile = candidate;
            return;
        end
    end

    candidates = struct([]);
    for iFolder = 1:numel(searchFolders)
        folderCandidates = dir(fullfile( ...
            searchFolders{iFolder}, 'part1_da3211_1_preprocessed*.mat'));
        if isempty(candidates)
            candidates = folderCandidates;
        else
            candidates = [candidates; folderCandidates]; %#ok<AGROW>
        end
    end

    if isempty(candidates)
        error([ ...
            'No Part 1 preprocessed MAT-file was found. Place it in the ' ...
            'same folder as this script or in the current MATLAB folder.']);
    end

    [~, newestIndex] = max([candidates.datenum]);
    dataFile = fullfile(candidates(newestIndex).folder, ...
                        candidates(newestIndex).name);

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

    figure('Name', 'Nominal representative errors with 3-sigma bounds');
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
        plot(t, 3*sigmaSignal, 'r--', 'LineWidth', 1.0);
        plot(t, -3*sigmaSignal, 'r--', 'LineWidth', 1.0);

        xlabel('Time [s]');
        ylabel(selectedLabels{iPlot});
        legend('Estimation error', '+3\sigma', '-3\sigma', ...
               'Location', 'best');
    end

    sgtitle('Nominal EKF representative errors and three-sigma bounds');

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


function plotNISComparison( ...
    t, NIS_nominal, NIS_highV, consistencyNominal, consistencyHighV)
%PLOTNISCOMPARISON Plot NIS and exact pointwise confidence limits.

    figure('Name', 'Normalized innovation squared comparison');
    set(gcf, 'Position', [100, 100, 1100, 520]);
    hold on;
    grid on;

    plot(t, NIS_nominal, 'r', 'LineWidth', 1.0);
    plot(t, NIS_highV, 'k--', 'LineWidth', 1.0);

    yline(consistencyNominal.expectedNIS, 'b--', ...
        'Expected NIS mean');
    yline(consistencyNominal.NISLower, 'm:', ...
        'Pointwise lower bound');
    yline(consistencyNominal.NISUpper, 'm:', ...
        'Pointwise upper bound');

    xlabel('Time [s]');
    ylabel('NIS');
    legend('Nominal', 'High-V-noise', 'Expected mean', ...
           'Confidence limits', 'Location', 'best');
    title(sprintf([ ...
        'NIS: nominal ANIS %.2f, high-V ANIS %.2f'], ...
        consistencyNominal.ANIS, consistencyHighV.ANIS));

end


function plotNEESComparison( ...
    t, NEES_nominal, NEES_highV, consistencyNominal, consistencyHighV)
%PLOTNEESCOMPARISON Plot full-state NEES consistency.

    figure('Name', 'Normalized estimation error squared comparison');
    set(gcf, 'Position', [100, 100, 1100, 520]);
    hold on;
    grid on;

    plot(t, NEES_nominal, 'r', 'LineWidth', 1.0);
    plot(t, NEES_highV, 'k--', 'LineWidth', 1.0);

    yline(consistencyNominal.expectedNEES, 'b--', ...
        'Expected NEES mean');
    yline(consistencyNominal.NEESLower, 'm:', ...
        'Pointwise lower bound');
    yline(consistencyNominal.NEESUpper, 'm:', ...
        'Pointwise upper bound');

    xlabel('Time [s]');
    ylabel('NEES');
    legend('Nominal', 'High-V-noise', 'Expected mean', ...
           'Confidence limits', 'Location', 'best');
    title(sprintf([ ...
        'NEES: nominal ANEES %.2f, high-V ANEES %.2f'], ...
        consistencyNominal.ANEES, consistencyHighV.ANEES));

end


function [data, t, z_m, u_IMU, Q_imu, R, x_true, ...
          lambda_true, wind_true] = validateAndExtractPart1Data(data)
%VALIDATEANDEXTRACTPART1DATA Check dimensions and extract variables.

    requiredVariables = {'t', 'z_m', 'u_IMU', 'Q', 'R', 'x_true'};
    for iVariable = 1:numel(requiredVariables)
        variableName = requiredVariables{iVariable};
        if ~isfield(data, variableName)
            error('The Part 1 MAT-file is missing variable "%s".', ...
                  variableName);
        end
    end

    t = data.t(:);
    N = numel(t);

    z_m    = orientSamplesByRow(data.z_m,   N, 'z_m');
    u_IMU  = orientSamplesByRow(data.u_IMU, N, 'u_IMU');
    x_true = orientSamplesByRow(data.x_true, N, 'x_true');

    Q_imu = data.Q;
    R     = data.R;

    if size(z_m, 2) ~= 12
        error('z_m must have 12 measurement columns.');
    end
    if size(u_IMU, 2) ~= 6
        error('u_IMU must have 6 IMU-input columns.');
    end
    if size(x_true, 2) ~= 18
        error('x_true must have 18 state columns.');
    end
    if ~isequal(size(Q_imu), [6, 6])
        error('Q must be 6-by-6.');
    end
    if ~isequal(size(R), [12, 12])
        error('R must be 12-by-12.');
    end
    if any(~isfinite(t)) || any(diff(t) <= 0)
        error('t must be finite and strictly increasing.');
    end
    if any(~isfinite(z_m(:))) || any(~isfinite(u_IMU(:))) || ...
            any(~isfinite(x_true(:)))
        error('The Part 1 data contain NaN or Inf values.');
    end
    if norm(Q_imu - Q_imu.', 'fro') > 1e-10 || ...
            min(eig(makeSymmetric(Q_imu))) < -1e-12
        error('Q must be symmetric positive semidefinite.');
    end
    if norm(R - R.', 'fro') > 1e-10 || ...
            min(eig(makeSymmetric(R))) <= 0
        error('R must be symmetric positive definite.');
    end

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

end


function gravityCheck = estimateEffectiveGravity( ...
    t, x_true, u_IMU, lambda_true)
%ESTIMATEEFFECTIVEGRAVITY Diagnose gravity using the simulated truth.
%
% This is a validation diagnostic, not an EKF measurement update. It uses
% the known simulated biases and state histories to assess whether the g
% used by the process model is compatible with the supplied trajectory.

    correctedIMU = u_IMU - lambda_true(:).';

    u = x_true(:, 4);
    v = x_true(:, 5);
    w = x_true(:, 6);
    phi   = x_true(:, 7);
    theta = x_true(:, 8);

    A_x = correctedIMU(:, 1);
    A_y = correctedIMU(:, 2);
    A_z = correctedIMU(:, 3);
    p = correctedIMU(:, 4);
    q = correctedIMU(:, 5);
    r = correctedIMU(:, 6);

    u_dot = gradient(u, t);
    v_dot = gradient(v, t);
    w_dot = gradient(w, t);

    denominatorY = sin(phi).*cos(theta);
    denominatorZ = cos(phi).*cos(theta);

    gFromY = (v_dot - A_y - p.*w + r.*u) ./ denominatorY;
    gFromZ = (w_dot - A_z - q.*u + p.*v) ./ denominatorZ;

    validY = abs(denominatorY) > 0.15 & isfinite(gFromY) & ...
             abs(gFromY) < 20;
    validZ = abs(denominatorZ) > 0.50 & isfinite(gFromZ) & ...
             abs(gFromZ) < 20;

    gravityCheck.gFromYMedian = median(gFromY(validY));
    gravityCheck.gFromZMedian = median(gFromZ(validZ));
    gravityCheck.gFromYRobustStd = robustStandardDeviation(gFromY(validY));
    gravityCheck.gFromZRobustStd = robustStandardDeviation(gFromZ(validZ));
    gravityCheck.numberOfYPoints = nnz(validY);
    gravityCheck.numberOfZPoints = nnz(validZ);

end


function value = robustStandardDeviation(x)
%ROBUSTSTANDARDDEVIATION Median-absolute-deviation scale estimate.

    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
        return;
    end
    centre = median(x);
    value = 1.4826 * median(abs(x - centre));

end


function comparisonTable = buildStateComparisonTable( ...
    x_true, resultsNominal, resultsHighV, stateNames, stateUnits)
%BUILDSTATECOMPARISONTABLE Final values and RMSE for both cases.

    displayScale = ones(18, 1);
    displayScale(7:9)   = 180/pi;
    displayScale(13:15) = 180/pi;

    trueFinal = x_true(end, :).'.*displayScale;
    nominalFinal = resultsNominal.x_hat_cor(end, :).'.*displayScale;
    highVFinal   = resultsHighV.x_hat_cor(end, :).'.*displayScale;

    nominalFinalError = resultsNominal.final_error_cor.*displayScale;
    highVFinalError   = resultsHighV.final_error_cor.*displayScale;
    nominalRMSE = resultsNominal.RMSE_cor.*displayScale;
    highVRMSE   = resultsHighV.RMSE_cor.*displayScale;
    RMSEratio   = highVRMSE ./ max(nominalRMSE, eps);

    comparisonTable = table( ...
        stateNames, stateUnits, trueFinal, nominalFinal, highVFinal, ...
        nominalFinalError, highVFinalError, nominalRMSE, highVRMSE, ...
        RMSEratio, ...
        'VariableNames', {'State', 'Unit', 'TrueFinal', ...
        'NominalFinal', 'HighVFinal', 'NominalFinalError', ...
        'HighVFinalError', 'NominalRMSE', 'HighVRMSE', ...
        'HighV_to_Nominal_RMSE'});

end


function stateTable = buildStateConsistencyTable( ...
    t, results, stateNames, stateUnits, settings)
%BUILDSTATECONSISTENCYTABLE Quantify covariance consistency after burn-in.

    scale = ones(18, 1);
    scale(7:9)   = 180/pi;
    scale(13:15) = 180/pi;

    finalError = results.final_error_cor.*scale;
    finalSigma = results.sigma_x_cor(:, end).*scale;
    finalAbsZ  = abs(finalError) ./ max(finalSigma, eps);

    use = t >= t(1) + settings.burnInSeconds;
    normalizedError = results.normalized_state_error(use, :);

    coverage1 = mean(abs(normalizedError) <= 1, 1).';
    coverage2 = mean(abs(normalizedError) <= 2, 1).';
    coverage3 = mean(abs(normalizedError) <= 3, 1).';
    meanSquaredNormalizedError = mean(normalizedError.^2, 1, 'omitnan').';

    stateTable = table(stateNames, stateUnits, finalError, finalSigma, ...
        finalAbsZ, meanSquaredNormalizedError, coverage1, coverage2, ...
        coverage3, ...
        'VariableNames', {'State', 'Unit', 'FinalError', 'FinalSigma', ...
        'FinalAbsZ', 'MeanSquaredNormalizedError', 'Coverage1Sigma', ...
        'Coverage2Sigma', 'Coverage3Sigma'});

end


function nuisanceTable = buildNuisanceConsistencyTable( ...
    resultsNominal, resultsHighV, lambda_true, wind_true)
%BUILDNUISANCECONSISTENCYTABLE Validate biases and wind against truth.

    names = {'lambda_Ax'; 'lambda_Ay'; 'lambda_Az'; ...
             'lambda_p'; 'lambda_q'; 'lambda_r'; ...
             'W_x'; 'W_y'; 'W_z'};
    units = {'m/s^2'; 'm/s^2'; 'm/s^2'; ...
             'deg/s'; 'deg/s'; 'deg/s'; ...
             'm/s'; 'm/s'; 'm/s'};

    truth = nan(9, 1);
    if numel(lambda_true) == 6
        truth(1:6) = lambda_true;
    end
    if numel(wind_true) == 3
        truth(7:9) = wind_true;
    end

    scale = ones(9, 1);
    scale(4:6) = 180/pi;

    nominalFinal = resultsNominal.x_hat_cor(end, 10:18).'.*scale;
    highVFinal   = resultsHighV.x_hat_cor(end, 10:18).'.*scale;
    truthDisplay = truth.*scale;

    nominalError = nominalFinal - truthDisplay;
    highVError   = highVFinal - truthDisplay;
    nominalSigma = resultsNominal.sigma_x_cor(10:18, end).*scale;
    highVSigma   = resultsHighV.sigma_x_cor(10:18, end).*scale;

    nominalAbsZ = abs(nominalError) ./ max(nominalSigma, eps);
    highVAbsZ   = abs(highVError) ./ max(highVSigma, eps);

    nominalWithin3Sigma = nominalAbsZ <= 3;
    highVWithin3Sigma   = highVAbsZ <= 3;

    nuisanceTable = table(names, units, truthDisplay, ...
        nominalFinal, nominalError, nominalSigma, nominalAbsZ, ...
        nominalWithin3Sigma, highVFinal, highVError, highVSigma, ...
        highVAbsZ, highVWithin3Sigma, ...
        'VariableNames', {'State', 'Unit', 'TrueValue', ...
        'NominalFinal', 'NominalError', 'NominalSigma', ...
        'NominalAbsZ', 'NominalWithin3Sigma', ...
        'HighVFinal', 'HighVError', 'HighVSigma', ...
        'HighVAbsZ', 'HighVWithin3Sigma'});

end


function innovationTable = buildInnovationStatistics( ...
    t, normalizedInnovation, measurementNames, settings)
%BUILDINNOVATIONSTATISTICS Mean, variance, coverage and whiteness tests.

    use = t >= t(1) + settings.burnInSeconds;
    X = normalizedInnovation(:, use);
    nMeasurements = size(X, 1);

    sampleMean = zeros(nMeasurements, 1);
    sampleStd = zeros(nMeasurements, 1);
    lag1Correlation = zeros(nMeasurements, 1);
    fractionWithin3Sigma = zeros(nMeasurements, 1);
    ljungBoxStatistic = zeros(nMeasurements, 1);
    ljungBoxPValue = zeros(nMeasurements, 1);

    for iMeasurement = 1:nMeasurements
        x = X(iMeasurement, :).';
        x = x(isfinite(x));

        sampleMean(iMeasurement) = mean(x);
        sampleStd(iMeasurement) = std(x, 0);
        fractionWithin3Sigma(iMeasurement) = mean(abs(x) <= 3);

        if numel(x) >= 3
            C = corrcoef(x(1:end-1), x(2:end));
            lag1Correlation(iMeasurement) = C(1, 2);
        else
            lag1Correlation(iMeasurement) = NaN;
        end

        [ljungBoxStatistic(iMeasurement), ...
         ljungBoxPValue(iMeasurement)] = ...
            ljungBoxTestLocal(x, settings.maxInnovationLag);
    end

    whitenessAccepted = ljungBoxPValue >= ...
        (1 - settings.confidenceLevel);

    innovationTable = table(measurementNames, sampleMean, sampleStd, ...
        lag1Correlation, fractionWithin3Sigma, ljungBoxStatistic, ...
        ljungBoxPValue, whitenessAccepted, ...
        'VariableNames', {'Measurement', 'Mean', 'Std', ...
        'Lag1Correlation', 'FractionWithin3Sigma', ...
        'LjungBoxStatistic', 'LjungBoxPValue', ...
        'WhitenessAccepted'});

end


function consistency = analyseFilterConsistency(t, results, settings)
%ANALYSEFILTERCONSISTENCY NIS and NEES consistency after burn-in.

    use = t >= t(1) + settings.burnInSeconds;
    alpha = 1 - settings.confidenceLevel;
    nUsed = nnz(use);
    nZ = 12;
    nX = 18;

    consistency.expectedNIS = nZ;
    consistency.NISLower = chiSquareInverse(alpha/2, nZ);
    consistency.NISUpper = chiSquareInverse(1 - alpha/2, nZ);
    consistency.ANIS = mean(results.NIS(use), 'omitnan');
    consistency.ANISLower = ...
        chiSquareInverse(alpha/2, nUsed*nZ) / nUsed;
    consistency.ANISUpper = ...
        chiSquareInverse(1 - alpha/2, nUsed*nZ) / nUsed;
    consistency.NISFractionInside = mean( ...
        results.NIS(use) >= consistency.NISLower & ...
        results.NIS(use) <= consistency.NISUpper, 'omitnan');

    consistency.expectedNEES = nX;
    consistency.NEESLower = chiSquareInverse(alpha/2, nX);
    consistency.NEESUpper = chiSquareInverse(1 - alpha/2, nX);
    consistency.ANEES = mean(results.NEES(use), 'omitnan');
    consistency.ANEESLower = ...
        chiSquareInverse(alpha/2, nUsed*nX) / nUsed;
    consistency.ANEESUpper = ...
        chiSquareInverse(1 - alpha/2, nUsed*nX) / nUsed;
    consistency.NEESFractionInside = mean( ...
        results.NEES(use) >= consistency.NEESLower & ...
        results.NEES(use) <= consistency.NEESUpper, 'omitnan');

    consistency.numberOfSamples = nUsed;
    consistency.confidenceLevel = settings.confidenceLevel;

end


function printConsistencySummary(caseName, consistency)
%PRINTCONSISTENCYSUMMARY Print NIS and NEES report-ready values.

    fprintf('\n================ %s CONSISTENCY ================\n', ...
        caseName);
    fprintf([ ...
        'ANIS: %.4f; expected %.1f; aggregate %.1f%% interval ' ...
        '[%.4f, %.4f]\n'], ...
        consistency.ANIS, consistency.expectedNIS, ...
        100*consistency.confidenceLevel, ...
        consistency.ANISLower, consistency.ANISUpper);
    fprintf('Fraction of pointwise NIS values inside interval: %.2f %%\n', ...
        100*consistency.NISFractionInside);

    fprintf('Time-average NEES: %.4f; pointwise expected value %.1f\n', ...
        consistency.ANEES, consistency.expectedNEES);
    fprintf('Fraction of pointwise NEES values inside interval: %.2f %%\n', ...
        100*consistency.NEESFractionInside);
    fprintf([ ...
        'Important: this is one trajectory, not an ensemble of independent ' ...
        'Monte Carlo runs. Consecutive state errors are correlated, so the ' ...
        'very narrow aggregate ANEES interval is not used as a formal ' ...
        'acceptance test. Use the pointwise NEES, state-wise normalized ' ...
        'errors and 3-sigma coverage to diagnose inconsistency.\n']);

end


function q = chiSquareInverse(probability, degreesOfFreedom)
%CHISQUAREINVERSE Chi-square inverse without Statistics Toolbox.

    q = 2 * gammaincinv(probability, degreesOfFreedom/2);

end


function [Qstat, pValue] = ljungBoxTestLocal(x, maxLag)
%LJUNGBOXTESTLOCAL Ljung-Box innovation-whiteness test.

    x = x(isfinite(x));
    x = x - mean(x);
    N = numel(x);
    maxLag = min(maxLag, max(N - 2, 1));

    if N < 5 || sum(x.^2) <= eps
        Qstat = NaN;
        pValue = NaN;
        return;
    end

    rho = sampleAutocorrelation(x, maxLag);
    lag = (1:maxLag).';
    Qstat = N*(N + 2) * sum((rho(2:end).^2) ./ (N - lag));
    pValue = gammainc(Qstat/2, maxLag/2, 'upper');

end


function rho = sampleAutocorrelation(x, maxLag)
%SAMPLEAUTOCORRELATION Biased sample autocorrelation from lag 0.

    x = x(:) - mean(x);
    denominator = sum(x.^2);
    rho = zeros(maxLag + 1, 1);
    rho(1) = 1;

    if denominator <= eps
        return;
    end

    for lag = 1:maxLag
        rho(lag + 1) = ...
            sum(x(1:end-lag).*x(1+lag:end)) / denominator;
    end

end


function value = safeQuadraticForm(errorVector, covarianceMatrix)
%SAFEQUADRATICFORM Compute e' P^{-1} e robustly.

    covarianceMatrix = makeSymmetric(covarianceMatrix);
    if rcond(covarianceMatrix) > 1e-12
        value = errorVector.' * (covarianceMatrix \ errorVector);
    else
        value = errorVector.' * pinv(covarianceMatrix) * errorVector;
    end

end


function plotNuisanceErrorsWith3Sigma( ...
    t, results, lambda_true, wind_true, figureTitle)
%PLOTNUISANCEERRORSWITH3SIGMA Plot all 6 biases and all 3 wind errors.

    truth = nan(9, 1);
    if numel(lambda_true) == 6
        truth(1:6) = lambda_true;
    end
    if numel(wind_true) == 3
        truth(7:9) = wind_true;
    end

    labels = {'lambda_{A_x}', 'lambda_{A_y}', 'lambda_{A_z}', ...
              'lambda_p', 'lambda_q', 'lambda_r', ...
              'W_x', 'W_y', 'W_z'};
    units = {'m/s^2', 'm/s^2', 'm/s^2', ...
             'deg/s', 'deg/s', 'deg/s', ...
             'm/s', 'm/s', 'm/s'};

    figure('Name', figureTitle);
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact');

    for i = 1:9
        nexttile;
        hold on;
        grid on;

        stateIndex = i + 9;
        estimate = results.x_hat_cor(:, stateIndex);
        sigma = results.sigma_x_cor(stateIndex, :).';
        errorSignal = estimate - truth(i);

        if i >= 4 && i <= 6
            errorSignal = errorSignal*180/pi;
            sigma = sigma*180/pi;
        end

        plot(t, errorSignal, 'b', 'LineWidth', 1.0);
        plot(t, 3*sigma, 'r--', 'LineWidth', 1.0);
        plot(t, -3*sigma, 'r--', 'LineWidth', 1.0);
        yline(0, 'k:');

        xlabel('Time [s]');
        ylabel(sprintf('%s error [%s]', labels{i}, units{i}));
        legend('Error', '+3 sigma', '-3 sigma', 'Location', 'best');
    end

    sgtitle(figureTitle);

end


function plotInnovationAutocorrelations( ...
    t, normalizedInnovation, measurementNames, settings, figureTitle)
%PLOTINNOVATIONAUTOCORRELATIONS Graphical innovation-whiteness check.
%
% figureTitle is optional so this helper also remains compatible with calls
% from earlier versions of the script that supplied only four inputs.

    if nargin < 5 || isempty(figureTitle)
        figureTitle = 'Normalized-innovation autocorrelation';
    end

    use = t >= t(1) + settings.burnInSeconds;
    normalizedInnovation = normalizedInnovation(:, use);
    N = size(normalizedInnovation, 2);
    approximate95 = 1.96/sqrt(N);

    figure('Name', figureTitle);
    set(gcf, 'Position', [100, 100, 1300, 850]);
    tiledlayout(4, 3, 'TileSpacing', 'compact');

    for i = 1:12
        nexttile;
        hold on;
        grid on;

        x = normalizedInnovation(i, :).';
        x = x(isfinite(x));
        rho = sampleAutocorrelation(x, settings.maxInnovationLag);
        stem(0:settings.maxInnovationLag, rho, ...
             'filled', 'MarkerSize', 3);
        yline(approximate95, 'r:');
        yline(-approximate95, 'r:');
        xlabel('Lag [samples]');
        ylabel('Autocorrelation');
        title(measurementNames{i}, 'Interpreter', 'none');
    end

    sgtitle(figureTitle);

end


function plotRMSEChange(comparisonTable)
%PLOTRMSECHANGE Show sensitivity to sigma_V = 5 m/s.

    figure('Name', 'RMSE sensitivity to airspeed measurement noise');
    set(gcf, 'Position', [100, 100, 1200, 520]);
    stateCategory = categorical(comparisonTable.State);
    stateCategory = reordercats(stateCategory, comparisonTable.State);
    bar(stateCategory, comparisonTable.HighV_to_Nominal_RMSE);
    grid on;
    yline(1, 'k--', 'No change');
    xlabel('State');
    ylabel('RMSE ratio: high-V / nominal');
    title('Effect of increasing true-airspeed noise to 5 m/s');

end


function saveAllFigures(outputFolder, resolution)
%SAVEALLFIGURES Save every open figure as FIG and PNG.

    figures = findall(groot, 'Type', 'figure');
    figures = flipud(figures);

    for i = 1:numel(figures)
        fig = figures(i);
        figureName = get(fig, 'Name');
        if isempty(figureName)
            figureName = sprintf('Figure_%02d', i);
        end

        safeName = regexprep(figureName, '[^a-zA-Z0-9_-]', '_');
        safeName = regexprep(safeName, '_+', '_');

        savefig(fig, fullfile(outputFolder, [safeName, '.fig']));

        pngFile = fullfile(outputFolder, [safeName, '.png']);
        try
            exportgraphics(fig, pngFile, 'Resolution', resolution);
        catch
            saveas(fig, pngFile);
        end
    end

end
