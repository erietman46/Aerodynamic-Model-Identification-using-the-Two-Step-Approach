clc;
clear;
close all;

% PART 2 - Extended Kalman Filter for flight-path reconstruction
%
% State vector:
% x = [xE yE zE u v w phi theta psi ...
%      lambda_Ax lambda_Ay lambda_Az ...
%      lambda_p lambda_q lambda_r Wx Wy Wz]'
%
% IMU input:
% u_IMU = [Ax_m Ay_m Az_m p_m q_m r_m]'
%
% Measurement vector:
% z = [x_GPS y_GPS z_GPS VN_GPS VE_GPS VD_GPS ...
%      phi_GPS theta_GPS psi_GPS V alpha beta]'

rng(2026);
g = 9.803;

%% Load the pre-processed Part 1 data

scriptFolder = fileparts(mfilename('fullpath'));
dataFiles = dir(fullfile(scriptFolder, ...
    'part1_da3211_1_preprocessed*.mat'));

if isempty(dataFiles)
    error(['Put part1_da3211_1_preprocessed.mat in the same folder ' ...
           'as this script.']);
end

load(fullfile(dataFiles(1).folder, dataFiles(1).name));

t = t(:);
N = length(t);

%% Case 1: normal airspeed noise, sigma_V = 0.2 m/s

fprintf('Running the normal EKF case...\n');
normal = runEKF(t, z_m, u_IMU, Q, R, g);

%% Case 2: increased airspeed noise, sigma_V = 5 m/s

% Use the same noise realization, but increase its standard deviation.
V_true = sqrt(x_true(:,4).^2 + x_true(:,5).^2 + x_true(:,6).^2);
originalAirspeedNoise = z_m(:,10) - V_true;

z_high = z_m;
z_high(:,10) = V_true + (5/0.2)*originalAirspeedNoise;

R_high = R;
R_high(10,10) = 5^2;

fprintf('Running the EKF with sigma_V = 5 m/s...\n');
high = runEKF(t, z_high, u_IMU, Q, R_high, g);

%% Calculate estimation errors and RMSE

normal.error = normal.xCorrected - x_true;
high.error   = high.xCorrected - x_true;

normal.error(:,7:9) = wrapAngle(normal.error(:,7:9));
high.error(:,7:9)   = wrapAngle(high.error(:,7:9));

normal.RMSE = sqrt(mean(normal.error.^2,1));
high.RMSE   = sqrt(mean(high.error.^2,1));

rad2deg = 180/pi;

stateNames = {'xE','yE','zE','u','v','w','phi','theta','psi', ...
              'lambda_Ax','lambda_Ay','lambda_Az', ...
              'lambda_p','lambda_q','lambda_r','Wx','Wy','Wz'};

stateUnits = {'m','m','m','m/s','m/s','m/s','deg','deg','deg', ...
    'm/s^2','m/s^2','m/s^2','deg/s','deg/s','deg/s','m/s','m/s','m/s'};

RMSEscale = ones(1,18);
RMSEscale([7 8 9 13 14 15]) = rad2deg;
normalRMSEdisplay = normal.RMSE.*RMSEscale;
highRMSEdisplay = high.RMSE.*RMSEscale;

fprintf('\nState RMSE comparison\n');
fprintf('State        Unit          Normal case       sigma_V = 5 m/s\n');
for i = 1:18
    fprintf('%-10s   %-7s      %12.6g       %12.6g\n', ...
        stateNames{i}, stateUnits{i}, ...
        normalRMSEdisplay(i), highRMSEdisplay(i));
end

RMSETable = table(stateNames.', stateUnits.', ...
    normalRMSEdisplay.', highRMSEdisplay.', ...
    high.RMSE.'./normal.RMSE.', ...
    'VariableNames', {'State','Unit','NormalRMSE','HighNoiseRMSE', ...
    'RMSEratio'});

%% Print final bias and wind estimates

fprintf('\nFinal accelerometer biases [m/s^2]\n');
fprintf('                 Ax            Ay            Az\n');
fprintf('True:       %10.6f    %10.6f    %10.6f\n', lambda_true(1:3));
fprintf('Normal:     %10.6f    %10.6f    %10.6f\n', ...
    normal.xCorrected(end,10:12));
fprintf('High noise: %10.6f    %10.6f    %10.6f\n', ...
    high.xCorrected(end,10:12));

fprintf('\nFinal gyro biases [deg/s]\n');
fprintf('                  p             q             r\n');
fprintf('True:       %10.6f    %10.6f    %10.6f\n', ...
    lambda_true(4:6)*rad2deg);
fprintf('Normal:     %10.6f    %10.6f    %10.6f\n', ...
    normal.xCorrected(end,13:15)*rad2deg);
fprintf('High noise: %10.6f    %10.6f    %10.6f\n', ...
    high.xCorrected(end,13:15)*rad2deg);

fprintf('\nFinal wind estimates [m/s]\n');
fprintf('                 Wx            Wy            Wz\n');
fprintf('True:       %10.6f    %10.6f    %10.6f\n', wind_true);
fprintf('Normal:     %10.6f    %10.6f    %10.6f\n', ...
    normal.xCorrected(end,16:18));
fprintf('High noise: %10.6f    %10.6f    %10.6f\n', ...
    high.xCorrected(end,16:18));

%% Check whether the final bias and wind estimates are consistent

trueNuisance = [lambda_true; wind_true];
normalFinalNuisance = normal.xCorrected(end,10:18).';
highFinalNuisance = high.xCorrected(end,10:18).';

normalFinalSigma = normal.sigmaCorrected(end,10:18).';
highFinalSigma = high.sigmaCorrected(end,10:18).';

normalFinalError = normalFinalNuisance - trueNuisance;
highFinalError = highFinalNuisance - trueNuisance;

normalSigmaRatio = normalFinalError./normalFinalSigma;
highSigmaRatio = highFinalError./highFinalSigma;

nuisanceNamesTable = {'lambda_Ax';'lambda_Ay';'lambda_Az'; ...
    'lambda_p';'lambda_q';'lambda_r';'Wx';'Wy';'Wz'};

biasWindTable = table(nuisanceNamesTable, trueNuisance, ...
    normalFinalNuisance, highFinalNuisance, ...
    normalFinalError, highFinalError, ...
    normalFinalSigma, highFinalSigma, ...
    normalSigmaRatio, highSigmaRatio, ...
    'VariableNames', {'State','TrueValue','NormalEstimate', ...
    'HighNoiseEstimate','NormalError','HighNoiseError', ...
    'NormalSigma','HighNoiseSigma','NormalErrorOverSigma', ...
    'HighNoiseErrorOverSigma'});

fprintf('\nFinal bias and wind consistency table (SI units)\n');
disp(biasWindTable);

if any(abs(normalSigmaRatio) > 3)
    fprintf(['WARNING: The following nominal estimates are outside their ' ...
        '3-sigma bounds:\n']);
    disp(nuisanceNamesTable(abs(normalSigmaRatio) > 3));
end

fprintf('\nMean NIS (expected value = 12)\n');
fprintf('Normal case:     %.4f\n', mean(normal.NIS(2:end)));
fprintf('High-noise case: %.4f\n', mean(high.NIS(2:end)));

% For 12 measurements, these are the 95 percent chi-square limits.
NISlower = 4.4038;
NISupper = 23.3367;

% 95 percent limits for the average NIS with 6000 innovations.
ANISlower = 11.8764;
ANISupper = 12.1243;

normalNISfraction = mean(normal.NIS(2:end) >= NISlower & ...
                         normal.NIS(2:end) <= NISupper);
highNISfraction = mean(high.NIS(2:end) >= NISlower & ...
                       high.NIS(2:end) <= NISupper);

fprintf('Normal NIS fraction inside 95%% limits: %.2f%%\n', ...
    100*normalNISfraction);
fprintf('High-noise NIS fraction inside 95%% limits: %.2f%%\n', ...
    100*highNISfraction);

normalMeanNIS = mean(normal.NIS(2:end));
highMeanNIS = mean(high.NIS(2:end));

fprintf('95%% average-NIS limits: %.4f to %.4f\n', ...
    ANISlower, ANISupper);
fprintf('Normal average NIS is inside limits: %d\n', ...
    normalMeanNIS >= ANISlower && normalMeanNIS <= ANISupper);
fprintf('High-noise average NIS is inside limits: %d\n', ...
    highMeanNIS >= ANISlower && highMeanNIS <= ANISupper);

%% Numerical innovation statistics

innovationMean = mean(normal.normalizedInnovation(2:end,:),1);
innovationStandardDeviation = ...
    std(normal.normalizedInnovation(2:end,:),0,1);
innovationMaximumAutocorrelation = zeros(1,12);

for i = 1:12
    innovationMaximumAutocorrelation(i) = maxInnovationAutocorrelation( ...
        normal.normalizedInnovation(2:end,i),25);
end

innovationNamesTable = {'xGPS';'yGPS';'zGPS';'VN_GPS';'VE_GPS'; ...
    'VD_GPS';'phi_GPS';'theta_GPS';'psi_GPS';'V';'alpha';'beta'};

innovationTable = table(innovationNamesTable, innovationMean.', ...
    innovationStandardDeviation.', innovationMaximumAutocorrelation.', ...
    'VariableNames', {'Measurement','Mean','StandardDeviation', ...
    'MaximumAbsoluteAutocorrelation'});

fprintf('\nNormalized innovation statistics\n');
disp(innovationTable);

%% Figure 1: raw and filtered measurements

measurementNames = {'x GPS [m]','y GPS [m]','z GPS [m]', ...
    'V_N GPS [m/s]','V_E GPS [m/s]','V_D GPS [m/s]', ...
    'phi GPS [deg]','theta GPS [deg]','psi GPS [deg]', ...
    'V [m/s]','alpha [deg]','beta [deg]'};

angleMeasurements = [7 8 9 11 12];

figure('Name','Raw and filtered measurements');
for i = 1:12
    subplot(4,3,i);
    hold on;
    grid on;

    rawSignal = z_m(:,i);
    filteredSignal = normal.zCorrected(:,i);

    if ismember(i, angleMeasurements)
        rawSignal = rawSignal*rad2deg;
        filteredSignal = filteredSignal*rad2deg;
    end

    plot(t, rawSignal, 'k.', 'MarkerSize', 2);
    plot(t, filteredSignal, 'r', 'LineWidth', 1);
    title(measurementNames{i});
    xlabel('Time [s]');
end
legend('Raw measurement','Filtered measurement');
sgtitle('Normal case: raw and EKF-filtered measurements');

%% Figure 1b: difference between raw and filtered measurements

figure('Name','Raw minus filtered measurements');
for i = 1:12
    subplot(4,3,i);
    hold on;
    grid on;

    measurementDifference = z_m(:,i) - normal.zCorrected(:,i);

    if ismember(i, angleMeasurements)
        measurementDifference = wrapAngle(measurementDifference)*rad2deg;
    end

    plot(t, measurementDifference, 'b');
    yline(0, 'k:');
    title(measurementNames{i});
    xlabel('Time [s]');
end
sgtitle('Normal case: raw measurement minus filtered measurement');

%% Figure 2: true states and EKF estimates for both cases

mainStateNames = {'x_E [m]','y_E [m]','z_E [m]', ...
    'u [m/s]','v [m/s]','w [m/s]', ...
    'phi [deg]','theta [deg]','psi [deg]'};

figure('Name','Aircraft-state estimates');
for i = 1:9
    subplot(3,3,i);
    hold on;
    grid on;

    trueSignal = x_true(:,i);
    normalSignal = normal.xCorrected(:,i);
    highSignal = high.xCorrected(:,i);

    if i >= 7
        trueSignal = trueSignal*rad2deg;
        normalSignal = normalSignal*rad2deg;
        highSignal = highSignal*rad2deg;
    end

    plot(t, trueSignal, 'b', 'LineWidth', 1.2);
    plot(t, normalSignal, 'r--', 'LineWidth', 1);
    plot(t, highSignal, 'k:', 'LineWidth', 1);
    title(mainStateNames{i});
    xlabel('Time [s]');
end
legend('True','Normal EKF','sigma_V = 5 m/s');
sgtitle('True aircraft states and EKF estimates');

%% Figure 2b: direct effect of increased airspeed noise on all states

figure('Name','Effect of increased airspeed noise');
for i = 1:18
    subplot(6,3,i);
    hold on;
    grid on;

    estimateDifference = high.xCorrected(:,i) - normal.xCorrected(:,i);

    if i >= 7 && i <= 9
        estimateDifference = wrapAngle(estimateDifference)*rad2deg;
    elseif i >= 13 && i <= 15
        estimateDifference = estimateDifference*rad2deg;
    end

    plot(t, estimateDifference, 'b');
    yline(0, 'k:');
    title(stateNames{i});
    xlabel('Time [s]');
end
sgtitle(['State-estimate change caused by increasing ' ...
    'airspeed noise to 5 m/s']);

%% Figure 3: estimated IMU biases and wind

nuisanceNames = {'lambda Ax [m/s^2]','lambda Ay [m/s^2]', ...
    'lambda Az [m/s^2]','lambda p [deg/s]', ...
    'lambda q [deg/s]','lambda r [deg/s]', ...
    'W_x [m/s]','W_y [m/s]','W_z [m/s]'};

figure('Name','Bias and wind estimates');
for i = 1:9
    stateIndex = i + 9;
    subplot(3,3,i);
    hold on;
    grid on;

    normalSignal = normal.xCorrected(:,stateIndex);
    highSignal = high.xCorrected(:,stateIndex);
    trueValue = trueNuisance(i);

    if stateIndex >= 13 && stateIndex <= 15
        normalSignal = normalSignal*rad2deg;
        highSignal = highSignal*rad2deg;
        trueValue = trueValue*rad2deg;
    end

    plot(t, normalSignal, 'r', 'LineWidth', 1);
    plot(t, highSignal, 'k--', 'LineWidth', 1);
    yline(trueValue, 'b:', 'LineWidth', 1.2);
    title(nuisanceNames{i});
    xlabel('Time [s]');
end
legend('Normal EKF','sigma_V = 5 m/s','True value');
sgtitle('Estimated accelerometer biases, gyro biases and wind');

%% Figure 4: bias and wind errors with 3-sigma bounds

figure('Name','Bias and wind convergence');
for i = 1:9
    stateIndex = i + 9;
    subplot(3,3,i);
    hold on;
    grid on;

    errorSignal = normal.error(:,stateIndex);
    sigmaSignal = normal.sigmaCorrected(:,stateIndex);

    if stateIndex >= 13 && stateIndex <= 15
        errorSignal = errorSignal*rad2deg;
        sigmaSignal = sigmaSignal*rad2deg;
    end

    plot(t, errorSignal, 'b');
    plot(t, 3*sigmaSignal, 'r--');
    plot(t, -3*sigmaSignal, 'r--');
    title(nuisanceNames{i});
    xlabel('Time [s]');
end
legend('Estimation error','+3 sigma','-3 sigma');
sgtitle('Normal case: bias and wind convergence');

%% Figure 4b: zoomed bias and wind convergence

zoomStart = 2;
zoomIndex = t >= zoomStart;

figure('Name','Bias and wind convergence zoom');
for i = 1:9
    stateIndex = i + 9;
    subplot(3,3,i);
    hold on;
    grid on;

    errorSignal = normal.error(:,stateIndex);
    sigmaSignal = normal.sigmaCorrected(:,stateIndex);

    if stateIndex >= 13 && stateIndex <= 15
        errorSignal = errorSignal*rad2deg;
        sigmaSignal = sigmaSignal*rad2deg;
    end

    plot(t(zoomIndex), errorSignal(zoomIndex), 'b');
    plot(t(zoomIndex), 3*sigmaSignal(zoomIndex), 'r--');
    plot(t(zoomIndex), -3*sigmaSignal(zoomIndex), 'r--');
    title(nuisanceNames{i});
    xlabel('Time [s]');
end
legend('Estimation error','+3 sigma','-3 sigma');
sgtitle('Normal case: zoomed bias and wind convergence');

%% Figure 5: main-state convergence with 3-sigma bounds

figure('Name','State convergence');
for i = 1:9
    subplot(3,3,i);
    hold on;
    grid on;

    errorSignal = normal.error(:,i);
    sigmaSignal = normal.sigmaCorrected(:,i);

    if i >= 7
        errorSignal = errorSignal*rad2deg;
        sigmaSignal = sigmaSignal*rad2deg;
    end

    plot(t, errorSignal, 'b');
    plot(t, 3*sigmaSignal, 'r--');
    plot(t, -3*sigmaSignal, 'r--');
    title(mainStateNames{i});
    xlabel('Time [s]');
end
legend('Estimation error','+3 sigma','-3 sigma');
sgtitle('Normal case: state convergence');

%% Figure 5b: zoomed main-state convergence

figure('Name','State convergence zoom');
for i = 1:9
    subplot(3,3,i);
    hold on;
    grid on;

    errorSignal = normal.error(:,i);
    sigmaSignal = normal.sigmaCorrected(:,i);

    if i >= 7
        errorSignal = errorSignal*rad2deg;
        sigmaSignal = sigmaSignal*rad2deg;
    end

    plot(t(zoomIndex), errorSignal(zoomIndex), 'b');
    plot(t(zoomIndex), 3*sigmaSignal(zoomIndex), 'r--');
    plot(t(zoomIndex), -3*sigmaSignal(zoomIndex), 'r--');
    title(mainStateNames{i});
    xlabel('Time [s]');
end
legend('Estimation error','+3 sigma','-3 sigma');
sgtitle('Normal case: zoomed state convergence');

%% Figure 6: normalized innovations

figure('Name','Normalized innovations');
for i = 1:12
    subplot(4,3,i);
    hold on;
    grid on;

    plot(t(2:end), normal.normalizedInnovation(2:end,i), 'b');
    yline(3, 'r--');
    yline(-3, 'r--');
    yline(0, 'k:');
    title(measurementNames{i});
    xlabel('Time [s]');
end
sgtitle('Normal case: normalized innovations with 3-sigma bounds');

%% Figure 7: NIS convergence

figure('Name','NIS');
hold on;
grid on;
plot(t(2:end), normal.NIS(2:end), 'b');
yline(12, 'r--', 'Expected mean');
yline(normalMeanNIS, 'm-', 'Observed mean');
yline(NISlower, 'k--', '95% lower limit');
yline(NISupper, 'k--', '95% upper limit');
xlabel('Time [s]');
ylabel('NIS');
title('Normal case: normalized innovation squared');

%% Save the Part 2 results

save(fullfile(scriptFolder, 'part2_da3211_EKF_results.mat'), ...
    't', 'x_true', 'z_m', 'z_high', 'u_IMU', ...
    'Q', 'R', 'R_high', 'normal', 'high', ...
    'lambda_true', 'wind_true', 'RMSETable', 'biasWindTable', ...
    'innovationTable', 'NISlower', 'NISupper', ...
    'ANISlower', 'ANISupper', 'normalNISfraction', ...
    'highNISfraction', 'normalMeanNIS', 'highMeanNIS');

fprintf('\nPart 2 is complete. Results were saved to:\n');
fprintf('part2_da3211_EKF_results.mat\n');


%% =====================================================================
%  EXTENDED KALMAN FILTER
%  =====================================================================

function result = runEKF(t, z, u_IMU, Q, R, g)

    N = length(t);
    numberOfStates = 18;
    numberOfMeasurements = 12;

    % Initial state estimate
    x = zeros(numberOfStates,1);
    x(1:3) = z(1,1:3).';

    V0 = z(1,10);
    alpha0 = z(1,11);
    beta0 = z(1,12);

    x(4) = V0*cos(alpha0)*cos(beta0);
    x(5) = V0*sin(beta0);
    x(6) = V0*sin(alpha0)*cos(beta0);
    x(7:9) = z(1,7:9).';

    % Initial covariance
    degree = pi/180;

    initialStandardDeviation = [ ...
        10; 10; 10; ...
        5; 5; 5; ...
        2*degree; 2*degree; 2*degree; ...
        0.1; 0.1; 0.1; ...
        0.05*degree; 0.05*degree; 0.05*degree; ...
        20; 20; 20];

    P = diag(initialStandardDeviation.^2);

    % Storage
    result.xPredicted = zeros(N,numberOfStates);
    result.xCorrected = zeros(N,numberOfStates);
    result.PPredicted = zeros(numberOfStates,numberOfStates,N);
    result.PCorrected = zeros(numberOfStates,numberOfStates,N);
    result.sigmaCorrected = zeros(N,numberOfStates);
    result.zPredicted = zeros(N,numberOfMeasurements);
    result.zCorrected = zeros(N,numberOfMeasurements);
    result.innovation = zeros(N,numberOfMeasurements);
    result.normalizedInnovation = zeros(N,numberOfMeasurements);
    result.NIS = zeros(N,1);

    result.xPredicted(1,:) = x.';
    result.xCorrected(1,:) = x.';
    result.PPredicted(:,:,1) = P;
    result.PCorrected(:,:,1) = P;
    result.sigmaCorrected(1,:) = sqrt(diag(P)).';
    result.zPredicted(1,:) = measurementModel(x).';
    result.zCorrected(1,:) = measurementModel(x).';

    for k = 1:N-1

        dt = t(k+1) - t(k);
        imu = u_IMU(k,:).';

        % -------------------------------------------------------------
        % Step 1: one-step-ahead state prediction
        % -------------------------------------------------------------
        xPredicted = rk4Step(x, imu, dt, g);
        xPredicted(7:9) = wrapAngle(xPredicted(7:9));

        % -------------------------------------------------------------
        % Step 2: calculate the Jacobians F and H
        % -------------------------------------------------------------
        F = numericalJacobian( ...
            @(state) stateModel(state,imu,g), xPredicted);

        H = numericalJacobian(@measurementModel, xPredicted);

        % -------------------------------------------------------------
        % Step 3: discretize the linearized model
        % -------------------------------------------------------------
        G = noiseInputMatrix(xPredicted);
        [Phi,Gamma] = c2d(F,G,dt);

        % -------------------------------------------------------------
        % Step 4: covariance prediction
        % -------------------------------------------------------------
        PPredicted = Phi*P*Phi.' + Gamma*Q*Gamma.';
        PPredicted = (PPredicted + PPredicted.')/2;

        % Predicted measurement and innovation
        zPredicted = measurementModel(xPredicted);
        innovation = z(k+1,:).' - zPredicted;

        angleIndex = [7 8 9 11 12];
        innovation(angleIndex) = wrapAngle(innovation(angleIndex));

        innovationCovariance = H*PPredicted*H.' + R;
        innovationCovariance = ...
            (innovationCovariance + innovationCovariance.')/2;

        % -------------------------------------------------------------
        % Step 5: Kalman gain
        % -------------------------------------------------------------
        K = (PPredicted*H.')/innovationCovariance;

        % -------------------------------------------------------------
        % Step 6: measurement correction
        % -------------------------------------------------------------
        x = xPredicted + K*innovation;
        x(7:9) = wrapAngle(x(7:9));

        % -------------------------------------------------------------
        % Step 7: covariance correction
        % -------------------------------------------------------------
        I = eye(numberOfStates);
        P = (I-K*H)*PPredicted*(I-K*H).' + K*R*K.';
        P = (P + P.')/2;

        % Store the results
        result.xPredicted(k+1,:) = xPredicted.';
        result.xCorrected(k+1,:) = x.';
        result.PPredicted(:,:,k+1) = PPredicted;
        result.PCorrected(:,:,k+1) = P;
        result.sigmaCorrected(k+1,:) = ...
            sqrt(max(diag(P),0)).';
        result.zPredicted(k+1,:) = zPredicted.';
        result.zCorrected(k+1,:) = measurementModel(x).';
        result.innovation(k+1,:) = innovation.';
        result.normalizedInnovation(k+1,:) = ...
            (innovation./sqrt(max(diag(innovationCovariance),eps))).';
        result.NIS(k+1) = ...
            innovation.'*(innovationCovariance\innovation);
    end
end


%% Aircraft state equations

function xdot = stateModel(x, imu, g)

    u = x(4);
    v = x(5);
    w = x(6);

    phi = x(7);
    theta = x(8);
    psi = x(9);

    % Remove the estimated IMU biases from the measured IMU signals.
    Ax = imu(1) - x(10);
    Ay = imu(2) - x(11);
    Az = imu(3) - x(12);
    p  = imu(4) - x(13);
    q  = imu(5) - x(14);
    r  = imu(6) - x(15);

    Wx = x(16);
    Wy = x(17);
    Wz = x(18);

    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);

    xdot = zeros(18,1);

    % Position rates in the Earth frame
    xdot(1) = ...
        (u*ctheta + (v*sphi+w*cphi)*stheta)*cpsi ...
        - (v*cphi-w*sphi)*spsi + Wx;

    xdot(2) = ...
        (u*ctheta + (v*sphi+w*cphi)*stheta)*spsi ...
        + (v*cphi-w*sphi)*cpsi + Wy;

    xdot(3) = -u*stheta + (v*sphi+w*cphi)*ctheta + Wz;

    % Body-velocity rates
    xdot(4) = Ax - g*stheta + r*v - q*w;
    xdot(5) = Ay + g*ctheta*sphi + p*w - r*u;
    xdot(6) = Az + g*ctheta*cphi + q*u - p*v;

    % Euler-angle rates
    xdot(7) = p + q*sphi*tan(theta) + r*cphi*tan(theta);
    xdot(8) = q*cphi - r*sphi;
    xdot(9) = (q*sphi + r*cphi)/ctheta;

    % Biases and wind are modeled as constants.
    xdot(10:18) = 0;
end


%% Measurement equations

function z = measurementModel(x)

    u = x(4);
    v = x(5);
    w = x(6);

    phi = x(7);
    theta = x(8);
    psi = x(9);

    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);

    z = zeros(12,1);

    % GPS position
    z(1:3) = x(1:3);

    % GPS velocity in the Earth frame
    z(4) = ...
        (u*ctheta + (v*sphi+w*cphi)*stheta)*cpsi ...
        - (v*cphi-w*sphi)*spsi + x(16);

    z(5) = ...
        (u*ctheta + (v*sphi+w*cphi)*stheta)*spsi ...
        + (v*cphi-w*sphi)*cpsi + x(17);

    z(6) = -u*stheta + (v*sphi+w*cphi)*ctheta + x(18);

    % GPS attitude
    z(7:9) = x(7:9);

    % Airdata measurements
    z(10) = sqrt(u^2 + v^2 + w^2);
    z(11) = atan2(w,u);
    z(12) = atan2(v,sqrt(u^2+w^2));
end


%% One Runge-Kutta integration step

function xNext = rk4Step(x, imu, dt, g)

    k1 = stateModel(x,         imu, g);
    k2 = stateModel(x+dt*k1/2, imu, g);
    k3 = stateModel(x+dt*k2/2, imu, g);
    k4 = stateModel(x+dt*k3,   imu, g);

    xNext = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end


%% Numerical Jacobian

function J = numericalJacobian(fun, x)

    y = fun(x);
    J = zeros(length(y),length(x));

    for i = 1:length(x)
        step = 1e-6*max(1,abs(x(i)));

        xPlus = x;
        xMinus = x;
        xPlus(i) = xPlus(i) + step;
        xMinus(i) = xMinus(i) - step;

        J(:,i) = (fun(xPlus)-fun(xMinus))/(2*step);
    end
end


%% Matrix that maps IMU noise into the state equations

function G = noiseInputMatrix(x)

    u = x(4);
    v = x(5);
    w = x(6);
    phi = x(7);
    theta = x(8);

    G = zeros(18,6);

    % Accelerometer noise
    G(4,1) = 1;
    G(5,2) = 1;
    G(6,3) = 1;

    % p gyro noise
    G(5,4) =  w;
    G(6,4) = -v;
    G(7,4) =  1;

    % q gyro noise
    G(4,5) = -w;
    G(6,5) =  u;
    G(7,5) =  sin(phi)*tan(theta);
    G(8,5) =  cos(phi);
    G(9,5) =  sin(phi)/cos(theta);

    % r gyro noise
    G(4,6) =  v;
    G(5,6) = -u;
    G(7,6) =  cos(phi)*tan(theta);
    G(8,6) = -sin(phi);
    G(9,6) =  cos(phi)/cos(theta);
end


%% Wrap angles to the interval [-pi, pi)

function angle = wrapAngle(angle)
    angle = mod(angle+pi,2*pi)-pi;
end


%% Maximum absolute innovation autocorrelation for selected lags

function maximumValue = maxInnovationAutocorrelation(signal,maxLag)

    signal = signal - mean(signal);
    denominator = signal.'*signal;
    autocorrelation = zeros(maxLag,1);

    if denominator <= eps
        maximumValue = 0;
        return;
    end

    for lag = 1:maxLag
        autocorrelation(lag) = ...
            signal(1+lag:end).'*signal(1:end-lag)/denominator;
    end

    maximumValue = max(abs(autocorrelation));
end
