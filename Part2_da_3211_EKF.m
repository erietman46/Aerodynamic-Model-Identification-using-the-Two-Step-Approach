clc;
clear all;
close all;

% STATES:
% xE, yE, zE, u, v, w, phi, theta, psi,
% lambda_Ax, lambda_Ay, lambda_Az, lambda_p, lambda_q, lambda_r,
% Wx, Wy, Wz

% MEASURED OUTPUT:
% x_GPS, y_GPS, z_GPS, u_GPS, v_GPS, w_GPS,
% phi_GPS, theta_GPS, psi_GPS, V_m, alpha_m, beta_m

%_____________LOAD THE PREPROCESSED DATA__________________
dataFile = 'C:\Users\eltjo\OneDrive - Delft University of Technology\TU Delft\Specialization\Aerodynamic-Model-Identification-using-the-Two-Step-Approach\part1_da3211_1_preprocessed.mat';

% Fallback if the file is in the same folder as this script
if exist(dataFile, 'file') ~= 2
    dataFile = 'part1_da3211_1_preprocessed.mat';
end

preprocessed_variables = load(dataFile);

% Extract the variables from the loaded structure.
% This is necessary because preprocessed_variables = load(dataFile)
% stores the variables inside a structure.
t     = preprocessed_variables.t;
z_m   = preprocessed_variables.z_m;
u_IMU = preprocessed_variables.u_IMU;
Q     = preprocessed_variables.Q;
R     = preprocessed_variables.R;

% Extract true states and true parameters if they are available.
% These are used for validation of the EKF result.
if isfield(preprocessed_variables, 'x_true')
    x_true = preprocessed_variables.x_true;
else
    error('The preprocessed file does not contain x_true, so true-state comparison is impossible.');
end

if isfield(preprocessed_variables, 'lambda_true')
    lambda_true = preprocessed_variables.lambda_true;
end

if isfield(preprocessed_variables, 'wind_true')
    wind_true = preprocessed_variables.wind_true;
end

% Make sure the time vector is a column vector.
t = t(:);

% If x_true was accidentally stored as 18 x N instead of N x 18, transpose it.
if size(x_true,1) ~= length(t) && size(x_true,2) == length(t)
    x_true = x_true.';
end

%__________________INITIALIZATION__________________________

nx = 18;
nz = 12;
nu = 6;

N  = length(t);
dt = 0.01;
d2r = pi/180;

% For consistency with the EKF demo notation:
n  = nx;       % number of states
nm = nz;       % number of measurements
m  = nu;       % number of inputs

% The demo stores data as columns:
% X(:,k), Z(:,k), U(:,k)
Z_k = z_m.';
U_k = u_IMU.';

% Initial position states
x0 = ones(nx, 1);
x0(1) = z_m(1,1);
x0(2) = z_m(1,2);
x0(3) = z_m(1,3);

% Initial body-axis velocity states.
% These are reconstructed from V, alpha, beta.
V0     = z_m(1,10);
alpha0 = z_m(1,11);
beta0  = z_m(1,12);

x0(4) = V0*cos(alpha0)*cos(beta0);
x0(5) = V0*sin(beta0);
x0(6) = V0*sin(alpha0)*cos(beta0);

% Initial attitude angles
x0(7) = z_m(1,7);
x0(8) = z_m(1,8);
x0(9) = z_m(1,9);

% Initial IMU bias guesses
x0(10:15) = 0.00001;

% Initial wind guesses
x0(16:18) = 0.00001;

% Initial covariance matrix of the process noise: Q
Q_imu = Q;

% Initial covariance matrix of the measurement noise: R
R_meas = R;

% Initial covariance matrix of the prediction: P0,0
std0 = [...
    10; 10; 10; ...                       % position uncertainty
    5; 5; 5; ...                       % body velocity uncertainty
    2*d2r; 2*d2r; 2*d2r; ...           % attitude uncertainty
    0.1; 0.1; 0.1; ...                 % accelerometer bias uncertainty
    0.05*d2r; 0.05*d2r; 0.05*d2r; ...  % angular rates bias uncertainty
    20; 20; 20];                       % wind uncertainty

P0 = diag(std0.^2);

% Same notation as the EKF demo.
E_x_0 = x0;
P_0   = P0;

%_____________________RUN THE EKF______________________________

% Initialize the EKF
t_k             = 0;
t_k1            = dt;

% Corrected EKF state estimate x(k+1|k+1)
XX_k1_k1        = zeros(n, N);

% Predicted EKF state estimate x(k+1|k)
% This is added so that final EKF prediction can also be compared with x_true.
XX_k1_k         = zeros(n, N);

% Corrected covariance matrix P(k+1|k+1)
PP_k1_k1        = zeros(n, n, N);

% Predicted covariance matrix P(k+1|k)
PP_k1_k         = zeros(n, n, N);

% Corrected state standard deviation
STD_x_cor       = zeros(n, N);

% Predicted state standard deviation
STD_x_pred      = zeros(n, N);

% Predicted measurement standard deviation
STD_z           = zeros(nm, N);

% Predicted measurement vector
ZZ_pred         = zeros(nm, N);

% Innovation vector
Innov_k         = zeros(nm, N);

% Normalized innovation squared
NIS_k           = zeros(1, N);

x_k1_k1         = E_x_0;    % x(0|0) = E{x_0}
P_k1_k1         = P_0;      % P(0|0) = P(0)

% Store initial values
XX_k1_k1(:,1)      = x_k1_k1;
XX_k1_k(:,1)       = x_k1_k1;
PP_k1_k1(:,:,1)    = P_k1_k1;
PP_k1_k(:,:,1)     = P_k1_k1;
STD_x_cor(:,1)     = sqrt(diag(P_k1_k1));
STD_x_pred(:,1)    = sqrt(diag(P_k1_k1));
ZZ_pred(:,1)       = measurementModel12_local(x_k1_k1);
Innov_k(:,1)       = innovation12_local(Z_k(:,1), ZZ_pred(:,1));

tic;

% Run the filter through all N samples
for k = 2:N

    % Use the actual timestep from the time vector.
    dt = t(k) - t(k-1);
    t_k = t(k-1);
    t_k1 = t(k);

    % x(k+1|k), prediction
    [~, x_k1_k] = rk4_local(@navDynamics18_local, x_k1_k1, U_k(:,k-1), [t_k, t_k1]);

    % Wrap attitude angles after propagation.
    x_k1_k(7:9) = wrapPi_local(x_k1_k(7:9));

    % Calculate Jacobians Phi(k+1,k) and Gamma(k+1,k)

    % Fx is the perturbation matrix of f(x,u,t)
    Fx = numericalJacobian_local(@(x) navDynamics18_local(0, x, U_k(:,k-1)), x_k1_k);

    % G is the continuous input-noise distribution matrix.
    % It maps IMU noise into the state derivatives.
    G = imuNoiseMapping_local(x_k1_k, U_k(:,k-1));

    % Continuous-to-discrete transformation
    [Phi, Gamma] = c2d_local(Fx, G, dt);

    % P(k+1|k), prediction covariance matrix
    P_k1_k = Phi*P_k1_k1*Phi.' + Gamma*Q_imu*Gamma.';

    % Numerical symmetry protection
    P_k1_k = 0.5*(P_k1_k + P_k1_k.');

    % Store predicted state and covariance
    XX_k1_k(:,k)    = x_k1_k;
    PP_k1_k(:,:,k)  = P_k1_k;
    STD_x_pred(:,k) = sqrt(diag(P_k1_k));

    % Correction

    % Hx is the perturbation matrix of h(x,u,t)
    Hx = numericalJacobian_local(@measurementModel12_local, x_k1_k);

    % Observation and observation error predictions
    z_k1_k = measurementModel12_local(x_k1_k);

    % Innovation: measured output minus predicted output
    innov = innovation12_local(Z_k(:,k), z_k1_k);

    % Covariance matrix of observation error
    P_zz = Hx*P_k1_k*Hx.' + R_meas;

    % Standard deviation of observation error for validation
    std_z = sqrt(diag(P_zz));

    % K(k+1), Kalman gain
    K = P_k1_k * Hx.' / P_zz;

    % Calculate optimal state x(k+1|k+1)
    x_k1_k1 = x_k1_k + K*innov;

    % Wrap attitude angles after correction.
    x_k1_k1(7:9) = wrapPi_local(x_k1_k1(7:9));

    % P(k+1|k+1), correction
    % Numerically stable Joseph form
    P_k1_k1 = (eye(n) - K*Hx)*P_k1_k*(eye(n) - K*Hx).' + K*R_meas*K.';

    % Numerical symmetry protection
    P_k1_k1 = 0.5*(P_k1_k1 + P_k1_k1.');

    std_x_cor = sqrt(diag(P_k1_k1));

    % Store corrected results
    XX_k1_k1(:,k)      = x_k1_k1;
    PP_k1_k1(:,:,k)    = P_k1_k1;
    STD_x_cor(:,k)     = std_x_cor;
    STD_z(:,k)         = std_z;
    ZZ_pred(:,k)       = z_k1_k;
    Innov_k(:,k)       = innov;

    % Normalized innovation squared.
    % For a 12-dimensional measurement, the expected mean is around 12.
    NIS_k(k) = innov.' * (P_zz \ innov);

    if mod(k,1000) == 0 || k == N
        fprintf('EKF running: k = %d / %d, t = %.2f s\n', k, N, t(k));
    end

end

time = toc;

fprintf('\nEKF completed run with %d samples in %2.2f seconds.\n', N, time);

% Convert estimates to convenient N x 18 matrices.
xhat_pred = XX_k1_k.';
xhat      = XX_k1_k1.';

%% _____________________TRUE-STATE COMPARISON______________________________

% Estimation errors over the complete maneuver
EstErr_x_pred = xhat_pred - x_true;
EstErr_x_cor  = xhat      - x_true;

% Wrap attitude errors, because phi/theta/psi are angular states.
EstErr_x_pred(:,7:9) = wrapPi_local(EstErr_x_pred(:,7:9));
EstErr_x_cor(:,7:9)  = wrapPi_local(EstErr_x_cor(:,7:9));

% Final true, predicted, and corrected states
x_true_final = x_true(end,:).';
x_pred_final = xhat_pred(end,:).';
x_cor_final  = xhat(end,:).';

% Final errors
final_pred_error = x_pred_final - x_true_final;
final_cor_error  = x_cor_final  - x_true_final;

final_pred_error(7:9) = wrapPi_local(final_pred_error(7:9));
final_cor_error(7:9)  = wrapPi_local(final_cor_error(7:9));

% RMSE over the complete maneuver
RMSE_pred = sqrt(mean(EstErr_x_pred.^2, 1)).';
RMSE_cor  = sqrt(mean(EstErr_x_cor.^2, 1)).';

% Names, units, and scale factors for readable printing.
stateNames = { ...
    'x_E'; ...
    'y_E'; ...
    'z_E'; ...
    'u'; ...
    'v'; ...
    'w'; ...
    'phi'; ...
    'theta'; ...
    'psi'; ...
    'lambda_Ax'; ...
    'lambda_Ay'; ...
    'lambda_Az'; ...
    'lambda_p'; ...
    'lambda_q'; ...
    'lambda_r'; ...
    'W_x'; ...
    'W_y'; ...
    'W_z'};

stateUnits = { ...
    'm'; ...
    'm'; ...
    'm'; ...
    'm/s'; ...
    'm/s'; ...
    'm/s'; ...
    'deg'; ...
    'deg'; ...
    'deg'; ...
    'm/s^2'; ...
    'm/s^2'; ...
    'm/s^2'; ...
    'deg/s'; ...
    'deg/s'; ...
    'deg/s'; ...
    'm/s'; ...
    'm/s'; ...
    'm/s'};

scale = ones(18,1);
scale(7:9)   = 180/pi;   % Euler angles
scale(13:15) = 180/pi;   % gyro biases

TrueFinal      = x_true_final    .* scale;
EKFPredFinal   = x_pred_final    .* scale;
EKFCorrFinal   = x_cor_final     .* scale;
PredFinalError = final_pred_error .* scale;
CorrFinalError = final_cor_error  .* scale;
PredRMSE       = RMSE_pred       .* scale;
CorrRMSE       = RMSE_cor        .* scale;

comparisonTable = table( ...
    stateNames, ...
    stateUnits, ...
    TrueFinal, ...
    EKFPredFinal, ...
    EKFCorrFinal, ...
    PredFinalError, ...
    CorrFinalError, ...
    PredRMSE, ...
    CorrRMSE, ...
    'VariableNames', { ...
        'State', ...
        'Unit', ...
        'TrueFinal', ...
        'EKFPredictedFinal', ...
        'EKFCorrectedFinal', ...
        'PredictedFinalError', ...
        'CorrectedFinalError', ...
        'PredictedRMSE', ...
        'CorrectedRMSE'});

fprintf('\n================ TRUE-STATE COMPARISON TABLE ================\n');
disp(comparisonTable);

fprintf('\nImportant interpretation:\n');
fprintf('EKFPredictedFinal is x(k+1|k), before using the final measurement.\n');
fprintf('EKFCorrectedFinal is x(k+1|k+1), after using the final measurement.\n');
fprintf('For reporting, the corrected estimate is usually the main EKF result.\n');

fprintf('\n================ FINAL ESTIMATED BIASES AND WIND ================\n');

fprintf('\nFinal estimated accelerometer biases [m/s^2]:\n');
fprintf('lambda_Ax = %.8f\n', xhat(end,10));
fprintf('lambda_Ay = %.8f\n', xhat(end,11));
fprintf('lambda_Az = %.8f\n', xhat(end,12));

fprintf('\nFinal estimated gyro biases [rad/s]:\n');
fprintf('lambda_p = %.10f\n', xhat(end,13));
fprintf('lambda_q = %.10f\n', xhat(end,14));
fprintf('lambda_r = %.10f\n', xhat(end,15));

fprintf('\nFinal estimated gyro biases [deg/s]:\n');
fprintf('lambda_p = %.8f\n', xhat(end,13)*180/pi);
fprintf('lambda_q = %.8f\n', xhat(end,14)*180/pi);
fprintf('lambda_r = %.8f\n', xhat(end,15)*180/pi);

fprintf('\nFinal estimated wind components [m/s]:\n');
fprintf('Wx = %.8f\n', xhat(end,16));
fprintf('Wy = %.8f\n', xhat(end,17));
fprintf('Wz = %.8f\n', xhat(end,18));

fprintf('\nMean NIS after first second = %.4f\n', mean(NIS_k(t > 1), 'omitnan'));
fprintf('Expected mean NIS is approximately nz = 12.\n');

if exist('lambda_true','var')
    fprintf('\nTrue IMU biases:\n');
    disp(lambda_true(:).');

    fprintf('Bias estimation error, corrected final estimate:\n');
    disp((xhat(end,10:15).' - lambda_true(:)).');
end

if exist('wind_true','var')
    fprintf('\nTrue wind:\n');
    disp(wind_true(:).');

    fprintf('Wind estimation error, corrected final estimate:\n');
    disp((xhat(end,16:18).' - wind_true(:)).');
end

%% _____________________SAVE RESULTS______________________________

save('part2_da3211_EKF_results_extended.mat', ...
    't', ...
    'x_true', ...
    'xhat', ...
    'xhat_pred', ...
    'EstErr_x_cor', ...
    'EstErr_x_pred', ...
    'RMSE_cor', ...
    'RMSE_pred', ...
    'comparisonTable', ...
    'XX_k1_k1', ...
    'XX_k1_k', ...
    'PP_k1_k1', ...
    'PP_k1_k', ...
    'STD_x_cor', ...
    'STD_x_pred', ...
    'STD_z', ...
    'ZZ_pred', ...
    'Innov_k', ...
    'NIS_k', ...
    'Q_imu', ...
    'R_meas');

fprintf('\nSaved results to part2_da3211_EKF_results_extended.mat\n');

%% _____________________PLOTTING______________________________

% Measurement names
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

%% Plot raw measurements and predicted/filtered measurements

plotID = 1001;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(4,3,'TileSpacing','compact');

for i = 1:12

    nexttile;
    hold on;
    grid on;

    measuredSignal = Z_k(i,:);
    predictedSignal = ZZ_pred(i,:);

    if ismember(i, angleMeasurementIndex)
        measuredSignal = measuredSignal*180/pi;
        predictedSignal = predictedSignal*180/pi;
    end

    plot(t, measuredSignal, 'k.');
    plot(t, predictedSignal, 'r', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(measurementNames{i});
    legend('Measurement', 'EKF predicted measurement', 'Location', 'best');

end

sgtitle('Raw measurements and EKF predicted measurements');

%% Plot true state and corrected estimated state

plotID = 1002;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(3,3,'TileSpacing','compact');

stateNames_9 = { ...
    'x_E [m]', ...
    'y_E [m]', ...
    'z_E [m]', ...
    'u [m/s]', ...
    'v [m/s]', ...
    'w [m/s]', ...
    '\phi [deg]', ...
    '\theta [deg]', ...
    '\psi [deg]'};

for i = 1:9

    nexttile;
    hold on;
    grid on;

    trueSignal = x_true(:,i);
    estimatedSignal = xhat(:,i);

    if i >= 7
        trueSignal = trueSignal*180/pi;
        estimatedSignal = estimatedSignal*180/pi;
    end

    plot(t, trueSignal, 'b', 'LineWidth', 1.2);
    plot(t, estimatedSignal, 'r--', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(stateNames_9{i});
    legend('True state', 'EKF corrected estimate', 'Location', 'best');

end

sgtitle('True state and EKF corrected state');

%% Plot true state, predicted state, and corrected state

plotID = 1004;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(3,3,'TileSpacing','compact');

for i = 1:9

    nexttile;
    hold on;
    grid on;

    trueSignal = x_true(:,i);
    predictedSignal = xhat_pred(:,i);
    correctedSignal = xhat(:,i);

    if i >= 7
        trueSignal = trueSignal*180/pi;
        predictedSignal = predictedSignal*180/pi;
        correctedSignal = correctedSignal*180/pi;
    end

    plot(t, trueSignal, 'b', 'LineWidth', 1.2);
    plot(t, predictedSignal, 'k:', 'LineWidth', 1.0);
    plot(t, correctedSignal, 'r--', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(stateNames_9{i});
    legend('True', 'EKF prediction x(k+1|k)', 'EKF correction x(k+1|k+1)', 'Location', 'best');

end

sgtitle('True state versus EKF prediction and correction');

%% Plot estimated bias and wind states

plotID = 1003;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(3,3,'TileSpacing','compact');

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

for j = 1:9

    nexttile;
    hold on;
    grid on;

    estimatedSignal = xhat(:,estimateIndex(j));

    if j >= 4 && j <= 6
        estimatedSignal = estimatedSignal*180/pi;
    end

    plot(t, estimatedSignal, 'r', 'LineWidth', 1.2);

    if exist('lambda_true','var') && j <= 6

        trueValue = lambda_true(j);

        if j >= 4
            trueValue = trueValue*180/pi;
        end

        yline(trueValue, 'b--', 'LineWidth', 1.2);

        legend('Estimated', 'True', 'Location', 'best');

    elseif exist('wind_true','var') && j >= 7

        trueValue = wind_true(j-6);

        yline(trueValue, 'b--', 'LineWidth', 1.2);

        legend('Estimated', 'True', 'Location', 'best');

    else

        legend('Estimated', 'Location', 'best');

    end

    xlabel('Time [s]');
    ylabel(estimateNames{j});

end

sgtitle('Estimated IMU biases and wind states');

%% Plot state standard deviations

plotID = 2001;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(3,3,'TileSpacing','compact');

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

for j = 1:length(selectedStates)

    nexttile;
    hold on;
    grid on;

    stateIndex = selectedStates(j);

    stdSignal = STD_x_cor(stateIndex,:);

    if stateIndex == 7
        stdSignal = stdSignal*180/pi;
    end

    plot(t, stdSignal, 'r', 'LineWidth', 1.2);

    xlabel('Time [s]');
    ylabel(selectedNames{j});

end

sgtitle('State estimation standard deviations');

%% Plot corrected state estimation errors with +/- standard deviation

plotID = 2002;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(3,3,'TileSpacing','compact');

for j = 1:length(selectedStates)

    nexttile;
    hold on;
    grid on;

    stateIndex = selectedStates(j);

    errorSignal = EstErr_x_cor(:,stateIndex);
    stdSignal = STD_x_cor(stateIndex,:).';

    if stateIndex == 7
        errorSignal = errorSignal*180/pi;
        stdSignal = stdSignal*180/pi;
    end

    plot(t, errorSignal, 'b', 'LineWidth', 1.0);
    plot(t, stdSignal, 'r--', 'LineWidth', 1.0);
    plot(t, -stdSignal, 'r--', 'LineWidth', 1.0);

    xlabel('Time [s]');
    ylabel(selectedNames{j});
    legend('Estimation error', '+1\sigma', '-1\sigma', 'Location', 'best');

end

sgtitle('Corrected EKF state estimation error with standard deviation bounds');

%% Plot normalized innovations

plotID = 3001;
figure(plotID);
set(plotID, 'Position', [100 100 1200 800], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

tiledlayout(4,3,'TileSpacing','compact');

Rdiag = diag(R_meas);

for i = 1:12

    nexttile;
    hold on;
    grid on;

    normalizedInnovation = Innov_k(i,:) ./ sqrt(Rdiag(i));

    plot(t, normalizedInnovation, 'b', 'LineWidth', 1.0);
    yline(0, 'k-');
    yline(3, 'r:');
    yline(-3, 'r:');

    xlabel('Time [s]');
    ylabel(['\nu_', num2str(i), '/\sigma']);

end

sgtitle('Normalized innovation sequence');

%% Plot NIS consistency

plotID = 3002;
figure(plotID);
set(plotID, 'Position', [100 100 900 500], ...
    'defaultaxesfontsize', 10, ...
    'defaulttextfontsize', 10, ...
    'PaperPositionMode', 'auto');

hold on;
grid on;

plot(t, NIS_k, 'b', 'LineWidth', 1.0);

yline(12, 'k--', 'Expected mean, n_z = 12');
yline(21.03, 'r:', 'Approx. 95% chi-square bound');
yline(26.22, 'r:', 'Approx. 99% chi-square bound');

xlabel('Time [s]');
ylabel('NIS');
title('Normalized Innovation Squared');

%% ======================= LOCAL FUNCTIONS ============================
% MATLAB allows local functions at the end of a script.

function [tout, xout] = rk4_local(fhandle, x0, u, tspan)

    t0 = tspan(1);
    t1 = tspan(2);
    dt = t1 - t0;

    k1 = fhandle(t0,          x0,              u);
    k2 = fhandle(t0 + dt/2,   x0 + dt*k1/2,    u);
    k3 = fhandle(t0 + dt/2,   x0 + dt*k2/2,    u);
    k4 = fhandle(t1,          x0 + dt*k3,      u);

    xout = x0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    tout = t1;

end

function xdot = navDynamics18_local(~, x, u_m)

    g = 9.80665;

    % States
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

    % Measured IMU inputs
    Ax_m = u_m(1);
    Ay_m = u_m(2);
    Az_m = u_m(3);

    p_m = u_m(4);
    q_m = u_m(5);
    r_m = u_m(6);

    % Bias-corrected IMU inputs
    Ax = Ax_m - lambda_Ax;
    Ay = Ay_m - lambda_Ay;
    Az = Az_m - lambda_Az;

    p = p_m - lambda_p;
    q = q_m - lambda_q;
    r = r_m - lambda_r;

    % Trigonometric terms
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
            - (v*cphi - w*sphi)*spsi + Wx;

    xdot(2) = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
            + (v*cphi - w*sphi)*cpsi + Wy;

    xdot(3) = -u*stheta + (v*sphi + w*cphi)*ctheta + Wz;

    % Body velocity dynamics
    xdot(4) = Ax - g*stheta + r*v - q*w;
    xdot(5) = Ay + g*sphi*ctheta + p*w - r*u;
    xdot(6) = Az + g*cphi*ctheta + q*u - p*v;

    % Attitude dynamics
    xdot(7) = p + q*sphi*ttheta + r*cphi*ttheta;
    xdot(8) = q*cphi - r*sphi;
    xdot(9) = (q*sphi + r*cphi)/ctheta;

    % Bias and wind dynamics
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

    % GPS ground-speed components
    zhat(4) = (u*ctheta + (v*sphi + w*cphi)*stheta)*cpsi ...
            - (v*cphi - w*sphi)*spsi + Wx;

    zhat(5) = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
            + (v*cphi - w*sphi)*cpsi + Wy;

    zhat(6) = -u*stheta + (v*sphi + w*cphi)*ctheta + Wz;

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

function G = imuNoiseMapping_local(x, ~)

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

    % Accelerometer noise
    G(4,1) = 1;
    G(5,2) = 1;
    G(6,3) = 1;

    % Gyro p noise
    G(5,4) = w;
    G(6,4) = -v;
    G(7,4) = 1;

    % Gyro q noise
    G(4,5) = -w;
    G(6,5) = u;
    G(7,5) = sphi*ttheta;
    G(8,5) = cphi;
    G(9,5) = sphi/ctheta;

    % Gyro r noise
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

    Phi = B(1:n,1:n);
    Gamma = B(1:n,n+1:n+m);

end

function nu = innovation12_local(z, zhat)

    nu = z - zhat;

    % Wrap angle residuals:
    % phi, theta, psi, alpha, beta
    angleIndex = [7 8 9 11 12];

    nu(angleIndex) = wrapPi_local(nu(angleIndex));

end

function angle = wrapPi_local(angle)

    angle = mod(angle + pi, 2*pi) - pi;

end