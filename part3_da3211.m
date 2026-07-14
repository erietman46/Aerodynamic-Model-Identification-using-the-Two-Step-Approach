clear; clc; close all;

% PART 3 - AERODYNAMIC MODEL IDENTIFICATION USING OLS AND RLS
% Manoeuvre: da3211_1 (aileron doublet)
%
% This self-contained script addresses all six Part 3 requirements for the
% lateral-directional derivatives that are identifiable in da3211_1.
%
% DATA FILE ROLES
%   da3211_1*.mat:
%       independent raw manoeuvre file; supplies the commanded controls.
%   part2_da3211_EKF_results*.mat:
%       output of Step 1; supplies the EKF state trajectory, corrected IMU
%       inputs and estimated accelerometer/gyro biases.
% The raw simulated states and aerodynamic coefficients are deliberately
% not used as regression outputs. This keeps Part 3 a genuine two-step
% identification using the estimates obtained in Part 2.
%
% Linear-in-the-parameters model:
%   y = A*theta + epsilon
%
% Initial ordinary least-squares estimate:
%   theta_OLS = inv(A_0'*A_0)*A_0'*y_0
%
% Recursive update for each new regression row a_(k+1):
%   e_(k+1)     = y_(k+1)-a_(k+1)*theta_k
%   K_(k+1)     = P_k*a_(k+1)'/(lambda+a_(k+1)*P_k*a_(k+1)')
%   theta_(k+1) = theta_k+K_(k+1)*e_(k+1)
%   P_(k+1)     = (I-K_(k+1)*a_(k+1))*P_k/lambda
%
% lambda = 1 is used because the simulated aerodynamic derivatives are
% constant. With lambda = 1, final RLS must equal batch OLS applied to the
% initial and new data together (apart from numerical round-off).
%
% IMPORTANT MANOEUVRE-SPECIFIC MODEL CHOICE
% Only the aileron is excited. Elevator, rudder and throttle are constant.
% Consequently this record can identify the lateral-directional equations
%
%   CY = CY0 + CY_beta*beta + CY_p*p*b/(2V)
%            + CY_r*r*b/(2V) + CY_da*delta_a
%   Cl = Cl0 + Cl_beta*beta + Cl_p*p*b/(2V)
%            + Cl_r*r*b/(2V) + Cl_da*delta_a
%   Cn = Cn0 + Cn_beta*beta + Cn_p*p*b/(2V)
%            + Cn_r*r*b/(2V) + Cn_da*delta_a,
%
% but it cannot identify longitudinal derivatives or delta_r derivatives.
% Only CY, Cl and Cn are therefore estimated. All six force and moment
% coefficients are nevertheless reconstructed as a check on Step 1.
%
% DATA ALLOCATION (strictly chronological and non-overlapping)
%   Initial OLS: 0.5 s trim + first aileron pulse.
%   RLS update:  second pulse of the first doublet.
%   Validation:  second doublet + five seconds of free decay.
% Thus the final RLS estimate uses the complete first doublet, while the
% second doublet is never used for parameter estimation.

%% Settings

showFigures = false;         % Set true to display figures interactively.
derivativeWindow = 51;      % 0.51 s at the 100 Hz sample rate.
derivativeOrder = 3;        % Cubic local-polynomial differentiation.
forgettingFactor = 1.0;     % Constant aircraft parameters: no forgetting.
validationDecayTime = 5.0;  % Free decay retained after validation pulse.
maximumResidualLag = 50;    % Residual ACF up to 0.50 s.
hacLag = 50;                % Newey-West/HAC covariance lag (0.50 s).

%% Load original flight data and Part 2 EKF results

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

rawFile = findDataFile(scriptFolder,'da3211_1*.mat');
ekfFile = findDataFile(scriptFolder,'part2_da3211_EKF_results*.mat');

raw = load(rawFile);
ekf = load(ekfFile,'t','normal','u_IMU');

t = ekf.t(:);
N = numel(t);

if numel(raw.t) ~= N || size(ekf.normal.xCorrected,1) ~= N || ...
        size(ekf.u_IMU,1) ~= N
    error('Original data and EKF results have different record lengths.');
end
if any(abs(raw.t(:)-t) > 1e-10)
    error('Original data and EKF time vectors do not match.');
end

outputFolder = fullfile(scriptFolder,'part3_output');
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% Appendix-A aircraft and atmosphere data

aircraft.mass = 5000;       % kg
aircraft.Ixx = 15351;       % kg m^2
aircraft.Iyy = 22965;       % kg m^2
aircraft.Izz = 36220;       % kg m^2
aircraft.Ixz = 1908;        % kg m^2
aircraft.b = 13.3250;       % m
aircraft.S = 24.9900;       % m^2
aircraft.c = 1.9910;        % m

testAltitude = 1000;        % m, stated in da3211_1 description
rho = isaDensity(testAltitude);

%% Extract EKF states and bias-correct the IMU measurements

xHat = ekf.normal.xCorrected;

u = xHat(:,4);
v = xHat(:,5);
w = xHat(:,6);
V = sqrt(u.^2+v.^2+w.^2);
alpha = atan2(w,u);
beta = atan2(v,sqrt(u.^2+w.^2));

accelerometerBias = xHat(:,10:12);
gyroBias = xHat(:,13:15);

specificForce = ekf.u_IMU(:,1:3)-accelerometerBias;
bodyRates = ekf.u_IMU(:,4:6)-gyroBias;

p = bodyRates(:,1);
q = bodyRates(:,2);
r = bodyRates(:,3);

%% Numerical differentiation required for aerodynamic moments

% A direct two-point finite difference strongly amplifies rate-gyro noise.
% Instead, at each sample a cubic polynomial is least-squares fitted over a
% 51-sample window. The linear coefficient of this local polynomial is the
% derivative at the central sample. This is a Savitzky-Golay-type numerical
% differentiator, implemented below without requiring a toolbox.

pDot = localPolynomialDerivative(t,p,derivativeOrder,derivativeWindow);
qDot = localPolynomialDerivative(t,q,derivativeOrder,derivativeWindow);
rDot = localPolynomialDerivative(t,r,derivativeOrder,derivativeWindow);

%% Reconstruct aerodynamic forces, moments and coefficients

dynamicPressure = 0.5*rho*V.^2;

Xforce = aircraft.mass*specificForce(:,1);
Yforce = aircraft.mass*specificForce(:,2);
Zforce = aircraft.mass*specificForce(:,3);

Lmoment = aircraft.Ixx*pDot ...
    + q.*r*(aircraft.Izz-aircraft.Iyy) ...
    - (p.*q+rDot)*aircraft.Ixz;

Mmoment = aircraft.Iyy*qDot ...
    + r.*p*(aircraft.Ixx-aircraft.Izz) ...
    + (p.^2-r.^2)*aircraft.Ixz;

Nmoment = aircraft.Izz*rDot ...
    + p.*q*(aircraft.Iyy-aircraft.Ixx) ...
    + (q.*r-pDot)*aircraft.Ixz;

CX = Xforce./(dynamicPressure*aircraft.S);
CY = Yforce./(dynamicPressure*aircraft.S);
CZ = Zforce./(dynamicPressure*aircraft.S);
Cl = Lmoment./(dynamicPressure*aircraft.S*aircraft.b);
Cm = Mmoment./(dynamicPressure*aircraft.S*aircraft.c);
Cn = Nmoment./(dynamicPressure*aircraft.S*aircraft.b);

coefficientData = [CX,CY,CZ,Cl,Cm,Cn];
allCoefficientNames = {'CX','CY','CZ','Cl','Cm','Cn'};

%% Complete Appendix-A regression structures

deltaA = raw.da(:);
deltaE = raw.de(:);
deltaR = raw.dr(:);
Tc = 0.5*(raw.Tc1(:)+raw.Tc2(:));

qHat = q*aircraft.c./V;
pHat = p*aircraft.b./(2*V);
rHat = r*aircraft.b./(2*V);

% Complete structures from Appendix A. These are saved even though this
% manoeuvre does not excite all their columns.
A_CX_complete = [ones(N,1),alpha,alpha.^2,qHat,deltaE,Tc];
A_CZ_Cm_complete = [ones(N,1),alpha,qHat,deltaE,Tc];
A_lateral_complete = [ones(N,1),beta,pHat,rHat,deltaA,deltaR];

completeNames.CX = {'1','alpha','alpha^2','qc/V','delta_e','Tc'};
completeNames.CZ_Cm = {'1','alpha','qc/V','delta_e','Tc'};
completeNames.lateral = ...
    {'1','beta','pb/(2V)','rb/(2V)','delta_a','delta_r'};

% delta_r is constant and collinear with the intercept; omit it.
A = A_lateral_complete(:,1:5);
Y = [CY,Cl,Cn];
parameterNames = completeNames.lateral(1:5);
responseNames = {'CY','Cl','Cn'};

%% Chronological and independent use of the aileron-doublet data
%
% da3211_1 contains four signed pulse segments: the first two form the
% identification doublet and the last two form the validation doublet.
% The initial OLS estimate uses only trim and pulse 1. RLS then processes
% pulse 2 sample by sample. Pulses 3-4 and their free decay remain unseen.

baselineSamples = t < min(t(end),t(1)+5);
deltaABaseline = median(deltaA(baselineSamples));
deltaAExcursion = deltaA-deltaABaseline;
excitationThreshold = 0.10*max(abs(deltaAExcursion));

pulseRuns = signedPulseRuns(deltaAExcursion,excitationThreshold);
if size(pulseRuns,1) < 4
    error(['Four signed aileron pulse segments were expected in da3211_1. ' ...
           'Inspect delta_a before changing the automatic data masks.']);
end

firstPulse = pulseRuns(1,1):pulseRuns(1,2);
updatePulse = pulseRuns(2,1):pulseRuns(2,2);
validationPulses = pulseRuns(3,1):pulseRuns(4,2);

preTrimMask = t >= max(t(1),t(firstPulse(1))-0.5) & ...
    t < t(firstPulse(1));
initialMask = preTrimMask;
initialMask(firstPulse) = true;

updateMask = false(N,1);
updateMask(updatePulse) = true;

validationPulseMask = false(N,1);
validationPulseMask(validationPulses) = true;
lastValidationPulse = pulseRuns(4,2);
validationDecayMask = t > t(lastValidationPulse) & ...
    t <= min(t(end),t(lastValidationPulse)+validationDecayTime);
validationMask = validationPulseMask | validationDecayMask;

if any(initialMask & updateMask) || any(initialMask & validationMask) || ...
        any(updateMask & validationMask)
    error('Identification, RLS-update and validation masks overlap.');
end
if max(t(initialMask)) >= min(t(updateMask)) || ...
        max(t(updateMask)) >= min(t(validationMask))
    error('The automatically selected data blocks are not chronological.');
end

A0 = A(initialMask,:);
Y0 = Y(initialMask,:);
Anew = A(updateMask,:);
Ynew = Y(updateMask,:);
Aval = A(validationMask,:);
Yval = Y(validationMask,:);

if rank(A0) < size(A0,2)
    error('The initial OLS regression matrix is rank deficient.');
end

combinedMask = initialMask | updateMask;
Acombined = A(combinedMask,:);
Ycombined = Y(combinedMask,:);

%% Initial OLS, recursive LS and batch-equivalence check

numberOfResponses = size(Y,2);
numberOfParameters = size(A,2);

OLS = repmat(emptyOLSStructure(),numberOfResponses,1);
RLS = repmat(emptyRLSStructure(),numberOfResponses,1);
BatchOLS = repmat(emptyOLSStructure(),numberOfResponses,1);

for outputIndex = 1:numberOfResponses
    % Initial course-note OLS estimate using the first data block.
    OLS(outputIndex) = fitOLS(A0,Y0(:,outputIndex), ...
        maximumResidualLag,hacLag);

    % Course-demo RLS initialization and sample-by-sample update.
    initialInverseInformation = pinv(A0.'*A0);
    RLS(outputIndex) = runRLS(Anew,Ynew(:,outputIndex), ...
        OLS(outputIndex).theta,initialInverseInformation, ...
        forgettingFactor,Acombined,Ycombined(:,outputIndex), ...
        maximumResidualLag,hacLag);

    % With lambda=1 this batch result must equal the final RLS estimate.
    BatchOLS(outputIndex) = fitOLS(Acombined,Ycombined(:,outputIndex), ...
        maximumResidualLag,hacLag);
    RLS(outputIndex).batchDifference = RLS(outputIndex).theta ...
        - BatchOLS(outputIndex).theta;
    RLS(outputIndex).maximumBatchDifference = ...
        max(abs(RLS(outputIndex).batchDifference));

    OLS(outputIndex).validation = modelMetrics( ...
        Yval(:,outputIndex),Aval*OLS(outputIndex).theta, ...
        numberOfParameters,maximumResidualLag,Aval(:,2:end));
    RLS(outputIndex).validation = modelMetrics( ...
        Yval(:,outputIndex),Aval*RLS(outputIndex).theta, ...
        numberOfParameters,maximumResidualLag,Aval(:,2:end));
end

%% Statistical quality and identifiability diagnostics

standardizedRegressors = standardizeColumns(Acombined(:,2:end));
normalizedConditionNumber = cond(standardizedRegressors);
regressorCorrelation = corrcoef(Acombined(:,2:end));
VIF = diag(pinv(regressorCorrelation));

% Parameter correlation matrices expose statistical coupling that cannot
% be seen from the diagonal parameter variances alone.
for outputIndex = 1:numberOfResponses
    RLS(outputIndex).parameterCorrelation = ...
        covarianceToCorrelation(RLS(outputIndex).hacCovariance);
end

%% Model-term dominance and alternative reduced structure

% Influence is the standard deviation of each term a_j*theta_j over all
% estimation/update samples. This compares actual contribution to output,
% not coefficient magnitude (regressors have different physical scales).
% The least influential non-intercept term is removed from each response.

termContribution = zeros(numberOfParameters,numberOfResponses);
termFraction = zeros(numberOfParameters,numberOfResponses);
leastInfluentialIndex = zeros(numberOfResponses,1);
Reduced = repmat(emptyOLSStructure(),numberOfResponses,1);
reducedColumns = cell(numberOfResponses,1);

for outputIndex = 1:numberOfResponses
    individualTerms = Acombined.*RLS(outputIndex).theta.';
    termContribution(1,outputIndex) = abs(RLS(outputIndex).theta(1));
    for parameterIndex = 2:numberOfParameters
        termContribution(parameterIndex,outputIndex) = ...
            std(individualTerms(:,parameterIndex));
    end

    dynamicContributionSum = sum(termContribution(2:end,outputIndex));
    termFraction(2:end,outputIndex) = ...
        termContribution(2:end,outputIndex)/dynamicContributionSum;

    [~,relativeIndex] = min(termContribution(2:end,outputIndex));
    leastInfluentialIndex(outputIndex) = relativeIndex+1;
    keepColumns = setdiff(1:numberOfParameters, ...
        leastInfluentialIndex(outputIndex),'stable');
    reducedColumns{outputIndex} = keepColumns;

    Reduced(outputIndex) = fitOLS(Acombined(:,keepColumns), ...
        Ycombined(:,outputIndex),maximumResidualLag,hacLag);
    Reduced(outputIndex).validation = modelMetrics( ...
        Yval(:,outputIndex),Aval(:,keepColumns)*Reduced(outputIndex).theta, ...
        numel(keepColumns),maximumResidualLag,Aval(:,keepColumns(2:end)));
end

%% Print assignment-ready results

fprintf('\nPART 3 - OLS AND RLS AERODYNAMIC IDENTIFICATION\n');
fprintf('Original flight file: %s\n',rawFile);
fprintf('EKF result file:      %s\n',ekfFile);
fprintf(['Raw file role: commanded inputs only. EKF file role: estimated ' ...
    'states, biases and IMU data.\n']);
fprintf('Air density at %.0f m: %.6f kg/m^3\n',testAltitude,rho);
fprintf('Initial OLS samples: %d (%.2f to %.2f s)\n', ...
    nnz(initialMask),min(t(initialMask)),max(t(initialMask)));
fprintf('New RLS samples:     %d (%.2f to %.2f s)\n', ...
    nnz(updateMask),min(t(updateMask)),max(t(updateMask)));
fprintf('Validation samples:  %d (%.2f to %.2f s)\n', ...
    nnz(validationMask),min(t(validationMask)),max(t(validationMask)));
fprintf('rank(A0) = %d of %d\n',rank(A0),numberOfParameters);
fprintf('Normalized condition number = %.3f\n',normalizedConditionNumber);
fprintf('std(delta_r) = %.3e rad: delta_r derivative is not identifiable.\n', ...
    std(deltaR));
fprintf(['HAC/Newey-West uncertainty uses %d lags; it is reported as the ' ...
    'primary uncertainty because equation-error residuals may be serially ' ...
    'correlated.\n'],hacLag);

for outputIndex = 1:numberOfResponses
    fprintf('\n%s PARAMETERS\n',responseNames{outputIndex});
    fprintf('%-12s %14s %14s %14s %14s %14s\n', ...
        'Term','Initial OLS','Final RLS','HAC std.err','CI95 low','CI95 high');
    for parameterIndex = 1:numberOfParameters
        fprintf('%-12s %+14.6e %+14.6e %14.6e %+14.6e %+14.6e\n', ...
            parameterNames{parameterIndex}, ...
            OLS(outputIndex).theta(parameterIndex), ...
            RLS(outputIndex).theta(parameterIndex), ...
            RLS(outputIndex).hacStandardError(parameterIndex), ...
            RLS(outputIndex).hacConfidenceInterval(parameterIndex,1), ...
            RLS(outputIndex).hacConfidenceInterval(parameterIndex,2));
    end

    fprintf('Maximum |RLS - combined batch OLS| = %.3e\n', ...
        RLS(outputIndex).maximumBatchDifference);
    fprintf('Initial OLS validation: R2 = %+.4f, RMSE = %.6e\n', ...
        OLS(outputIndex).validation.R2, ...
        OLS(outputIndex).validation.RMSE);
    fprintf('Final RLS validation:   R2 = %+.4f, RMSE = %.6e\n', ...
        RLS(outputIndex).validation.R2, ...
        RLS(outputIndex).validation.RMSE);
    fprintf(['RLS residual mean = %+.3e, DW = %.3f, ' ...
        'max |ACF| = %.3f (95%% limit %.3f)\n'], ...
        RLS(outputIndex).validation.meanResidual, ...
        RLS(outputIndex).validation.DurbinWatson, ...
        RLS(outputIndex).validation.maxAbsoluteAutocorrelation, ...
        RLS(outputIndex).validation.autocorrelation95Limit);
    fprintf(['Ljung-Box p = %.3e, Jarque-Bera p = %.3e, ' ...
        'max |corr(residual,regressor)| = %.3f\n'], ...
        RLS(outputIndex).validation.LjungBoxPValue, ...
        RLS(outputIndex).validation.JarqueBeraPValue, ...
        RLS(outputIndex).validation.maxAbsoluteRegressorResidualCorrelation);

    certainty = abs(RLS(outputIndex).theta(2:end)) ...
        ./max(RLS(outputIndex).hacStandardError(2:end),eps);
    [~,mostCertainRelative] = max(certainty);
    [~,leastCertainRelative] = min(certainty);
    fprintf('Most certain dynamic term:  %s (HAC |t| = %.2f)\n', ...
        parameterNames{mostCertainRelative+1},certainty(mostCertainRelative));
    fprintf('Least certain dynamic term: %s (HAC |t| = %.2f)\n', ...
        parameterNames{leastCertainRelative+1},certainty(leastCertainRelative));
    fprintf('Least influential term: %s (%.2f%% of dynamic contribution)\n', ...
        parameterNames{leastInfluentialIndex(outputIndex)}, ...
        100*termFraction(leastInfluentialIndex(outputIndex),outputIndex));
    fprintf('Reduced validation: R2 = %+.4f, RMSE = %.6e\n', ...
        Reduced(outputIndex).validation.R2, ...
        Reduced(outputIndex).validation.RMSE);
    if RLS(outputIndex).validation.LjungBoxPValue < 0.05
        fprintf(['WARNING: %s validation residuals are not white; model ' ...
            'structure and/or the stochastic error model remain incomplete.\n'], ...
            responseNames{outputIndex});
    end
end

fprintf('\nVIF values for the dynamic regressors\n');
for parameterIndex = 2:numberOfParameters
    fprintf('  %-10s %.3f\n',parameterNames{parameterIndex}, ...
        VIF(parameterIndex-1));
end

%% Figures

if showFigures
    figureVisibility = 'on';
else
    figureVisibility = 'off';
end

fig1 = figure('Visible',figureVisibility,'Name','Reconstructed coefficients');
subplot(4,1,1); hold on;
plot(t,deltaA,'k','LineWidth',0.9);
plot(t(initialMask),deltaA(initialMask),'.','Color',[0 0.45 0.74]);
plot(t(updateMask),deltaA(updateMask),'.','Color',[0.85 0.33 0.10]);
plot(t(validationMask),deltaA(validationMask),'.','Color',[0.47 0.67 0.19]);
grid on; ylabel('\delta_a [rad]');
title('Aileron input and reconstructed lateral coefficients');
legend('All data','Initial OLS','RLS update','Validation', ...
    'Location','best');
subplot(4,1,2); plot(t,CY,'b'); grid on; ylabel('C_Y');
subplot(4,1,3); plot(t,Cl,'b'); grid on; ylabel('C_l');
subplot(4,1,4); plot(t,Cn,'b'); grid on; ylabel('C_n'); xlabel('Time [s]');
print(fig1,fullfile(outputFolder,'part3_reconstructed_coefficients.png'), ...
    '-dpng','-r150');

figDerivative = figure('Visible',figureVisibility, ...
    'Name','Angular-rate differentiation');
rateMatrix = [p,q,r];
derivativeMatrix = [pDot,qDot,rDot];
rateLabels = {'p','q','r'};
derivativeLabels = {'dp/dt','dq/dt','dr/dt'};
for rateIndex = 1:3
    subplot(3,2,2*rateIndex-1);
    plot(t,rateMatrix(:,rateIndex),'b'); grid on;
    ylabel([rateLabels{rateIndex} ' [rad/s]']);
    if rateIndex == 1
        title('Bias-corrected angular rates');
    end
    subplot(3,2,2*rateIndex);
    plot(t,derivativeMatrix(:,rateIndex),'r'); grid on;
    ylabel([derivativeLabels{rateIndex} ' [rad/s^2]']);
    if rateIndex == 1
        title(sprintf('Local cubic derivative, %d samples', ...
            derivativeWindow));
    end
end
subplot(3,2,5); xlabel('Time [s]');
subplot(3,2,6); xlabel('Time [s]');
print(figDerivative,fullfile(outputFolder, ...
    'part3_angular_rate_differentiation.png'),'-dpng','-r150');

figForces = figure('Visible',figureVisibility,'Name','Forces and moments');
dimensionalOutputs = [Xforce,Yforce,Zforce,Lmoment,Mmoment,Nmoment];
dimensionalLabels = {'X [N]','Y [N]','Z [N]', ...
    'L [N m]','M [N m]','N [N m]'};
for outputIndex = 1:6
    subplot(2,3,outputIndex);
    plot(t,dimensionalOutputs(:,outputIndex),'b'); grid on;
    ylabel(dimensionalLabels{outputIndex});
    if outputIndex > 3; xlabel('Time [s]'); end
end
sgtitle('Aerodynamic forces and moments reconstructed from EKF outputs');
print(figForces,fullfile(outputFolder,'part3_forces_moments.png'), ...
    '-dpng','-r150');

fig2 = figure('Visible',figureVisibility,'Name','OLS and RLS validation');
for outputIndex = 1:numberOfResponses
    subplot(3,1,outputIndex); hold on; grid on;
    keepColumns = reducedColumns{outputIndex};
    plot(t(validationMask),Yval(:,outputIndex),'k','LineWidth',1.2);
    plot(t(validationMask),Aval*OLS(outputIndex).theta,'b--','LineWidth',1.1);
    plot(t(validationMask),Aval*RLS(outputIndex).theta,'r-.','LineWidth',1.1);
    plot(t(validationMask), ...
        Aval(:,keepColumns)*Reduced(outputIndex).theta,'m:','LineWidth',1.2);
    ylabel(responseNames{outputIndex});
    if outputIndex == 1
        title('Independent validation data');
        legend('Reconstructed','Initial OLS','Final RLS','Reduced', ...
            'Location','best');
    end
end
xlabel('Time [s]');
print(fig2,fullfile(outputFolder,'part3_OLS_RLS_validation.png'), ...
    '-dpng','-r150');

fig3 = figure('Visible',figureVisibility,'Name','RLS convergence');
for outputIndex = 1:numberOfResponses
    subplot(3,1,outputIndex); hold on; grid on;
    plot(t(updateMask),RLS(outputIndex).thetaHistory,'LineWidth',1.0);
    ylabel(responseNames{outputIndex});
    if outputIndex == 1
        title('RLS parameter convergence during new data');
        legend(parameterNames,'Location','best');
    end
end
xlabel('Time [s]');
print(fig3,fullfile(outputFolder,'part3_RLS_parameter_convergence.png'), ...
    '-dpng','-r150');

figComparison = figure('Visible',figureVisibility, ...
    'Name','OLS and RLS parameter comparison');
for outputIndex = 1:numberOfResponses
    subplot(3,2,2*outputIndex-1);
    bar([OLS(outputIndex).theta,RLS(outputIndex).theta]); grid on;
    set(gca,'XTick',1:numberOfParameters,'XTickLabel',parameterNames);
    ylabel(responseNames{outputIndex});
    if outputIndex == 1
        title('Parameter estimates');
        legend('Initial OLS','Final RLS','Location','best');
    end

    subplot(3,2,2*outputIndex);
    standardizedChange = (RLS(outputIndex).theta-OLS(outputIndex).theta) ...
        ./max(RLS(outputIndex).hacStandardError,eps);
    bar(standardizedChange); hold on; grid on;
    plot([0.5,numberOfParameters+0.5],[1.96,1.96],'r--');
    plot([0.5,numberOfParameters+0.5],[-1.96,-1.96],'r--');
    set(gca,'XTick',1:numberOfParameters,'XTickLabel',parameterNames);
    ylabel('\Delta\theta / HAC SE');
    if outputIndex == 1
        title('Change after the RLS data block');
    end
end
print(figComparison,fullfile(outputFolder, ...
    'part3_OLS_RLS_parameter_comparison.png'),'-dpng','-r150');

fig4 = figure('Visible',figureVisibility,'Name','RLS residual analysis');
for outputIndex = 1:numberOfResponses
    validationResidual = Yval(:,outputIndex)-Aval*RLS(outputIndex).theta;
    [acf,lags] = sampleAutocorrelation(validationResidual,maximumResidualLag);
    acfLimit = 1.96/sqrt(numel(validationResidual));

    subplot(3,3,3*outputIndex-2);
    plot(t(validationMask),validationResidual,'b'); grid on;
    ylabel([responseNames{outputIndex} ' residual']);
    if outputIndex == 1; title('Validation residuals'); end

    subplot(3,3,3*outputIndex-1);
    stem(lags(2:end),acf(2:end),'filled'); hold on; grid on;
    plot([1,lags(end)],[acfLimit,acfLimit],'r--');
    plot([1,lags(end)],[-acfLimit,-acfLimit],'r--');
    ylabel([responseNames{outputIndex} ' ACF']);
    if outputIndex == 1; title('Residual autocorrelation'); end

    subplot(3,3,3*outputIndex);
    sortedResidual = sort(validationResidual);
    probability = ((1:numel(sortedResidual)).'-0.5)/numel(sortedResidual);
    theoreticalNormal = sqrt(2)*erfinv(2*probability-1);
    normalizedResidual = (sortedResidual-mean(sortedResidual)) ...
        /max(std(sortedResidual),eps);
    plot(theoreticalNormal,normalizedResidual,'b.'); hold on; grid on;
    normalBounds = [min(theoreticalNormal),max(theoreticalNormal)];
    plot(normalBounds,normalBounds,'r--');
    ylabel([responseNames{outputIndex} ' quantiles']);
    if outputIndex == 1; title('Normal Q-Q check'); end
end
subplot(3,3,7); xlabel('Time [s]');
subplot(3,3,8); xlabel('Lag [samples]');
subplot(3,3,9); xlabel('Normal quantiles');
print(fig4,fullfile(outputFolder,'part3_RLS_residual_analysis.png'), ...
    '-dpng','-r150');

fig5 = figure('Visible',figureVisibility,'Name','Term dominance');
bar(100*termFraction(2:end,:).'); grid on;
set(gca,'XTick',1:3,'XTickLabel',responseNames);
xlabel('Coefficient model'); ylabel('Dynamic contribution [%]');
legend(parameterNames(2:end),'Location','best');
title('Relative model-term influence');
print(fig5,fullfile(outputFolder,'part3_term_dominance.png'), ...
    '-dpng','-r150');

%% Save numerical results and report tables

dataUse.initialMask = initialMask;
dataUse.updateMask = updateMask;
dataUse.validationMask = validationMask;
dataUse.initialExplanation = ...
    'Half-second trim and pulse 1 used for initial OLS.';
dataUse.updateExplanation = ...
    'Pulse 2 processed chronologically sample-by-sample by RLS.';
dataUse.validationExplanation = ...
    'Pulses 3-4 and free decay reserved for independent validation only.';
dataUse.rawFile = rawFile;
dataUse.rawFileRole = 'Commanded controls and manoeuvre timing.';
dataUse.ekfFile = ekfFile;
dataUse.ekfFileRole = ...
    'Estimated state trajectory, IMU inputs and estimated sensor biases.';

regression.complete.A_CX = A_CX_complete;
regression.complete.A_CZ_Cm = A_CZ_Cm_complete;
regression.complete.A_lateral = A_lateral_complete;
regression.complete.names = completeNames;
regression.identified.A = A;
regression.identified.parameterNames = parameterNames;
regression.identified.responseNames = responseNames;

save(fullfile(outputFolder,'part3_da3211_OLS_RLS_results.mat'), ...
    't','rho','aircraft','V','alpha','beta','specificForce','bodyRates', ...
    'pDot','qDot','rDot','Xforce','Yforce','Zforce','Lmoment','Mmoment', ...
    'Nmoment','coefficientData','allCoefficientNames','regression', ...
    'dataUse','OLS','RLS','BatchOLS','Reduced','termContribution', ...
    'termFraction','leastInfluentialIndex','reducedColumns', ...
    'normalizedConditionNumber','regressorCorrelation','VIF', ...
    'derivativeWindow','derivativeOrder','maximumResidualLag','hacLag');

writeParameterCSV(fullfile(outputFolder,'part3_OLS_RLS_parameters.csv'), ...
    responseNames,parameterNames,OLS,RLS,BatchOLS);
writeValidationCSV(fullfile(outputFolder,'part3_validation_metrics.csv'), ...
    responseNames,OLS,RLS,Reduced,parameterNames,leastInfluentialIndex);
writeDominanceCSV(fullfile(outputFolder,'part3_term_dominance.csv'), ...
    responseNames,parameterNames,termContribution,termFraction);
writeIdentifiabilityCSV(fullfile(outputFolder, ...
    'part3_identifiability_diagnostics.csv'),parameterNames, ...
    responseNames,RLS,normalizedConditionNumber,VIF,regressorCorrelation);
writeDataUseCSV(fullfile(outputFolder,'part3_data_use.csv'), ...
    rawFile,ekfFile,t,initialMask,updateMask,validationMask);

fprintf('\nPart 3 results saved in:\n%s\n',outputFolder);

%% =====================================================================
% Local functions
% =====================================================================

function filePath = findDataFile(scriptFolder,pattern)
    searchFolders = {scriptFolder,fullfile(scriptFolder,'upload'), ...
        pwd,fullfile(pwd,'upload')};
    matches = [];
    for folderIndex = 1:numel(searchFolders)
        if exist(searchFolders{folderIndex},'dir')
            candidate = dir(fullfile(searchFolders{folderIndex},pattern));
            if ~isempty(candidate)
                matches = candidate;
                break;
            end
        end
    end
    if isempty(matches)
        error('Could not find a file matching %s.',pattern);
    end
    [~,newestIndex] = max([matches.datenum]);
    filePath = fullfile(matches(newestIndex).folder,matches(newestIndex).name);
end

function rho = isaDensity(altitude)
    T0 = 288.15;
    p0 = 101325;
    lapseRate = 0.0065;
    gasConstant = 287.05287;
    gravity = 9.80665;
    T = T0-lapseRate*altitude;
    pressure = p0*(T/T0)^(gravity/(gasConstant*lapseRate));
    rho = pressure/(gasConstant*T);
end

function derivative = localPolynomialDerivative(t,signal,order,windowLength)
    t = t(:);
    signal = signal(:);
    sampleCount = numel(t);
    if mod(windowLength,2) == 0 || windowLength <= order
        error('Differentiation window must be odd and exceed the order.');
    end
    halfWindow = floor(windowLength/2);
    derivative = zeros(sampleCount,1);

    for sampleIndex = 1:sampleCount
        firstIndex = max(1,sampleIndex-halfWindow);
        lastIndex = min(sampleCount,sampleIndex+halfWindow);
        if firstIndex == 1
            lastIndex = min(sampleCount,windowLength);
        elseif lastIndex == sampleCount
            firstIndex = max(1,sampleCount-windowLength+1);
        end

        tau = t(firstIndex:lastIndex)-t(sampleIndex);
        polynomialMatrix = zeros(numel(tau),order+1);
        for powerIndex = 0:order
            polynomialMatrix(:,powerIndex+1) = tau.^powerIndex;
        end
        coefficients = polynomialMatrix\signal(firstIndex:lastIndex);
        derivative(sampleIndex) = coefficients(2);
    end
end

function runs = signedPulseRuns(inputExcursion,threshold)
    % Return [startIndex,endIndex,sign] for each above-threshold segment.
    inputState = zeros(numel(inputExcursion),1);
    inputState(inputExcursion > threshold) = 1;
    inputState(inputExcursion < -threshold) = -1;
    starts = [1;find(diff(inputState) ~= 0)+1];
    ends = [starts(2:end)-1;numel(inputState)];
    keep = inputState(starts) ~= 0;
    runs = [starts(keep),ends(keep),inputState(starts(keep))];
end

function fit = emptyOLSStructure()
    fit.theta = [];
    fit.inverseInformation = [];
    fit.residualVariance = [];
    fit.covariance = [];
    fit.variance = [];
    fit.standardError = [];
    fit.confidenceInterval = [];
    fit.tStatistic = [];
    fit.pValueNormalApproximation = [];
    fit.hacCovariance = [];
    fit.hacVariance = [];
    fit.hacStandardError = [];
    fit.hacConfidenceInterval = [];
    fit.hacTStatistic = [];
    fit.hacPValueNormalApproximation = [];
    fit.prediction = [];
    fit.residual = [];
    fit.metrics = struct();
    fit.validation = struct();
end

function fit = fitOLS(A,y,maximumResidualLag,hacLag)
    y = y(:);
    sampleCount = size(A,1);
    parameterCount = size(A,2);

    % Course-note estimator. pinv is used instead of explicit inv for
    % numerical robustness, while giving the same result at full rank.
    inverseInformation = pinv(A.'*A);
    theta = inverseInformation*A.'*y;
    prediction = A*theta;
    residual = y-prediction;
    residualVariance = (residual.'*residual)/(sampleCount-parameterCount);
    covariance = residualVariance*inverseInformation;
    standardError = sqrt(max(diag(covariance),0));
    tStatistic = theta./max(standardError,eps);
    hacCovariance = neweyWestCovariance(A,residual,hacLag);
    hacStandardError = sqrt(max(diag(hacCovariance),0));
    hacTStatistic = theta./max(hacStandardError,eps);

    fit.theta = theta;
    fit.inverseInformation = inverseInformation;
    fit.residualVariance = residualVariance;
    fit.covariance = covariance;
    fit.variance = diag(covariance);
    fit.standardError = standardError;
    fit.confidenceInterval = [theta-1.96*standardError, ...
        theta+1.96*standardError];
    fit.tStatistic = tStatistic;
    fit.pValueNormalApproximation = erfc(abs(tStatistic)/sqrt(2));
    fit.hacCovariance = hacCovariance;
    fit.hacVariance = diag(hacCovariance);
    fit.hacStandardError = hacStandardError;
    fit.hacConfidenceInterval = [theta-1.96*hacStandardError, ...
        theta+1.96*hacStandardError];
    fit.hacTStatistic = hacTStatistic;
    fit.hacPValueNormalApproximation = ...
        erfc(abs(hacTStatistic)/sqrt(2));
    fit.prediction = prediction;
    fit.residual = residual;
    fit.metrics = modelMetrics(y,prediction,parameterCount, ...
        maximumResidualLag,A(:,2:end));
    % Keep the returned field set identical to emptyOLSStructure so this
    % result can be assigned into a preallocated MATLAB structure array.
    fit.validation = struct();
end

function result = emptyRLSStructure()
    result.theta = [];
    result.inverseInformation = [];
    result.residualVariance = [];
    result.covariance = [];
    result.variance = [];
    result.standardError = [];
    result.confidenceInterval = [];
    result.tStatistic = [];
    result.pValueNormalApproximation = [];
    result.hacCovariance = [];
    result.hacVariance = [];
    result.hacStandardError = [];
    result.hacConfidenceInterval = [];
    result.hacTStatistic = [];
    result.hacPValueNormalApproximation = [];
    result.thetaHistory = [];
    result.inverseInformationHistory = [];
    result.innovation = [];
    result.gain = [];
    result.prediction = [];
    result.residual = [];
    result.metrics = struct();
    result.validation = struct();
    result.batchDifference = [];
    result.maximumBatchDifference = [];
    result.parameterCorrelation = [];
end

function result = runRLS(Anew,ynew,initialTheta,initialP,lambda, ...
        Acombined,ycombined,maximumResidualLag,hacLag)
    ynew = ynew(:);
    theta = initialTheta;
    P = initialP;
    parameterCount = numel(theta);
    newSampleCount = size(Anew,1);

    thetaHistory = zeros(newSampleCount,parameterCount);
    inverseInformationHistory = zeros(parameterCount,parameterCount,newSampleCount);
    innovation = zeros(newSampleCount,1);
    gain = zeros(newSampleCount,parameterCount);

    for sampleIndex = 1:newSampleCount
        a = Anew(sampleIndex,:);
        y = ynew(sampleIndex);

        aP = a*P;
        K = (P*a.')/(lambda+aP*a.');
        e = y-a*theta;
        theta = theta+K*e;
        P = (eye(parameterCount)-K*a)*P/lambda;
        P = (P+P.')/2;

        thetaHistory(sampleIndex,:) = theta.';
        inverseInformationHistory(:,:,sampleIndex) = P;
        innovation(sampleIndex) = e;
        gain(sampleIndex,:) = K.';
    end

    prediction = Acombined*theta;
    residual = ycombined-prediction;
    residualVariance = (residual.'*residual) ...
        /(size(Acombined,1)-parameterCount);
    covariance = residualVariance*P;
    standardError = sqrt(max(diag(covariance),0));
    tStatistic = theta./max(standardError,eps);
    hacCovariance = neweyWestCovariance(Acombined,residual,hacLag);
    hacStandardError = sqrt(max(diag(hacCovariance),0));
    hacTStatistic = theta./max(hacStandardError,eps);

    result.theta = theta;
    result.inverseInformation = P;
    result.residualVariance = residualVariance;
    result.covariance = covariance;
    result.variance = diag(covariance);
    result.standardError = standardError;
    result.confidenceInterval = [theta-1.96*standardError, ...
        theta+1.96*standardError];
    result.tStatistic = tStatistic;
    result.pValueNormalApproximation = erfc(abs(tStatistic)/sqrt(2));
    result.hacCovariance = hacCovariance;
    result.hacVariance = diag(hacCovariance);
    result.hacStandardError = hacStandardError;
    result.hacConfidenceInterval = [theta-1.96*hacStandardError, ...
        theta+1.96*hacStandardError];
    result.hacTStatistic = hacTStatistic;
    result.hacPValueNormalApproximation = ...
        erfc(abs(hacTStatistic)/sqrt(2));
    result.thetaHistory = thetaHistory;
    result.inverseInformationHistory = inverseInformationHistory;
    result.innovation = innovation;
    result.gain = gain;
    result.prediction = prediction;
    result.residual = residual;
    result.metrics = modelMetrics(ycombined,prediction,parameterCount, ...
        maximumResidualLag,Acombined(:,2:end));
    % Keep the returned field set identical to emptyRLSStructure so this
    % result can be assigned into a preallocated MATLAB structure array.
    result.validation = struct();
    result.batchDifference = [];
    result.maximumBatchDifference = [];
    result.parameterCorrelation = [];
end

function metrics = modelMetrics(observed,predicted,parameterCount,maxLag,regressors)
    observed = observed(:);
    predicted = predicted(:);
    residual = observed-predicted;
    sampleCount = numel(observed);
    SSE = residual.'*residual;
    centredOutput = observed-mean(observed);
    SST = centredOutput.'*centredOutput;

    if SST > eps
        R2 = 1-SSE/SST;
        adjustedR2 = 1-(1-R2)*(sampleCount-1) ...
            /max(sampleCount-parameterCount,1);
    else
        R2 = NaN;
        adjustedR2 = NaN;
    end

    if SSE > eps && sampleCount > 1
        durbinWatson = sum(diff(residual).^2)/SSE;
    else
        durbinWatson = NaN;
    end

    [acf,~] = sampleAutocorrelation(residual,min(maxLag,sampleCount-1));
    if numel(acf) > 1
        maxAbsoluteAutocorrelation = max(abs(acf(2:end)));
    else
        maxAbsoluteAutocorrelation = NaN;
    end

    acfWithoutZero = acf(2:end);
    acfLags = (1:numel(acfWithoutZero)).';
    if isempty(acfWithoutZero)
        ljungBoxQ = NaN;
        ljungBoxPValue = NaN;
    else
        ljungBoxQ = sampleCount*(sampleCount+2)*sum( ...
            acfWithoutZero.^2./(sampleCount-acfLags));
        ljungBoxPValue = gammainc(ljungBoxQ/2, ...
            numel(acfWithoutZero)/2,'upper');
    end

    centredResidual = residual-mean(residual);
    residualStd = sqrt(mean(centredResidual.^2));
    if residualStd > eps
        residualSkewness = mean((centredResidual/residualStd).^3);
        residualKurtosis = mean((centredResidual/residualStd).^4);
        jarqueBera = sampleCount/6*(residualSkewness^2 + ...
            0.25*(residualKurtosis-3)^2);
        % A chi-square distribution with two degrees of freedom has the
        % survival function exp(-x/2), avoiding a toolbox dependency.
        jarqueBeraPValue = exp(-jarqueBera/2);
    else
        residualSkewness = NaN;
        residualKurtosis = NaN;
        jarqueBera = NaN;
        jarqueBeraPValue = NaN;
    end

    maxRegressorResidualCorrelation = ...
        maximumLaggedCorrelation(residual,regressors,maxLag);

    metrics.meanResidual = mean(residual);
    metrics.standardDeviationResidual = std(residual);
    metrics.RMSE = sqrt(mean(residual.^2));
    metrics.R2 = R2;
    metrics.adjustedR2 = adjustedR2;
    metrics.DurbinWatson = durbinWatson;
    metrics.maxAbsoluteAutocorrelation = maxAbsoluteAutocorrelation;
    metrics.autocorrelation95Limit = 1.96/sqrt(sampleCount);
    metrics.LjungBoxQ = ljungBoxQ;
    metrics.LjungBoxPValue = ljungBoxPValue;
    metrics.residualSkewness = residualSkewness;
    metrics.residualKurtosis = residualKurtosis;
    metrics.JarqueBera = jarqueBera;
    metrics.JarqueBeraPValue = jarqueBeraPValue;
    metrics.maxAbsoluteRegressorResidualCorrelation = ...
        maxRegressorResidualCorrelation;
    metrics.AIC = sampleCount*log(max(SSE/sampleCount,eps)) ...
        + 2*parameterCount;
    metrics.BIC = sampleCount*log(max(SSE/sampleCount,eps)) ...
        + parameterCount*log(sampleCount);
end

function covariance = neweyWestCovariance(A,residual,maxLag)
    % Heteroscedasticity-and-autocorrelation-consistent covariance using
    % Bartlett weights. It remains valid when residuals are serially
    % correlated, unlike the classical sigma^2*(A'*A)^(-1) covariance.
    sampleCount = size(A,1);
    parameterCount = size(A,2);
    maxLag = min(maxLag,sampleCount-1);
    score = A.*residual;
    meat = score.'*score;
    for lag = 1:maxLag
        weight = 1-lag/(maxLag+1);
        lagProduct = score(1+lag:end,:).'*score(1:end-lag,:);
        meat = meat+weight*(lagProduct+lagProduct.');
    end
    bread = pinv(A.'*A);
    covariance = sampleCount/max(sampleCount-parameterCount,1) ...
        *bread*meat*bread;
    covariance = (covariance+covariance.')/2;
end

function maximumCorrelation = maximumLaggedCorrelation(x,regressors,maxLag)
    if isempty(regressors)
        maximumCorrelation = NaN;
        return;
    end
    x = x(:)-mean(x);
    maximumCorrelation = 0;
    for columnIndex = 1:size(regressors,2)
        z = regressors(:,columnIndex)-mean(regressors(:,columnIndex));
        for lag = -maxLag:maxLag
            if lag >= 0
                xa = x(1+lag:end);
                za = z(1:end-lag);
            else
                xa = x(1:end+lag);
                za = z(1-lag:end);
            end
            denominator = sqrt((xa.'*xa)*(za.'*za));
            if denominator > eps
                maximumCorrelation = max(maximumCorrelation, ...
                    abs(xa.'*za/denominator));
            end
        end
    end
end

function standardized = standardizeColumns(matrix)
    columnMean = mean(matrix,1);
    columnStd = std(matrix,0,1);
    columnStd(columnStd < eps) = 1;
    standardized = (matrix-columnMean)./columnStd;
end

function correlation = covarianceToCorrelation(covariance)
    scale = sqrt(max(diag(covariance),eps));
    correlation = covariance./(scale*scale.');
end

function [acf,lags] = sampleAutocorrelation(signal,maxLag)
    signal = signal(:)-mean(signal);
    denominator = signal.'*signal;
    maxLag = max(0,min(maxLag,numel(signal)-1));
    lags = (0:maxLag).';
    acf = zeros(maxLag+1,1);
    if denominator <= eps
        acf(1) = 1;
        return;
    end
    acf(1) = 1;
    for lag = 1:maxLag
        acf(lag+1) = signal(1+lag:end).'*signal(1:end-lag)/denominator;
    end
end

function writeParameterCSV(filePath,responseNames,parameterNames,OLS,RLS,BatchOLS)
    fileID = fopen(filePath,'w');
    if fileID < 0; error('Could not create %s.',filePath); end
    fprintf(fileID,['Response,Parameter,Initial_OLS,Final_RLS,' ...
        'Combined_batch_OLS,RLS_variance,RLS_standard_error,' ...
        'RLS_CI95_low,RLS_CI95_high,RLS_t_statistic,RLS_p_normal,' ...
        'HAC_variance,HAC_standard_error,HAC_CI95_low,HAC_CI95_high,' ...
        'HAC_t_statistic,HAC_p_normal\n']);
    for outputIndex = 1:numel(responseNames)
        for parameterIndex = 1:numel(parameterNames)
            fprintf(fileID,['%s,%s,%.12g,%.12g,%.12g,%.12g,%.12g,' ...
                '%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,' ...
                '%.12g,%.12g\n'], ...
                responseNames{outputIndex},parameterNames{parameterIndex}, ...
                OLS(outputIndex).theta(parameterIndex), ...
                RLS(outputIndex).theta(parameterIndex), ...
                BatchOLS(outputIndex).theta(parameterIndex), ...
                RLS(outputIndex).variance(parameterIndex), ...
                RLS(outputIndex).standardError(parameterIndex), ...
                RLS(outputIndex).confidenceInterval(parameterIndex,1), ...
                RLS(outputIndex).confidenceInterval(parameterIndex,2), ...
                RLS(outputIndex).tStatistic(parameterIndex), ...
                RLS(outputIndex).pValueNormalApproximation(parameterIndex), ...
                RLS(outputIndex).hacVariance(parameterIndex), ...
                RLS(outputIndex).hacStandardError(parameterIndex), ...
                RLS(outputIndex).hacConfidenceInterval(parameterIndex,1), ...
                RLS(outputIndex).hacConfidenceInterval(parameterIndex,2), ...
                RLS(outputIndex).hacTStatistic(parameterIndex), ...
                RLS(outputIndex).hacPValueNormalApproximation(parameterIndex));
        end
    end
    fclose(fileID);
end

function writeValidationCSV(filePath,responseNames,OLS,RLS,Reduced, ...
        parameterNames,leastInfluentialIndex)
    fileID = fopen(filePath,'w');
    if fileID < 0; error('Could not create %s.',filePath); end
    fprintf(fileID,['Response,Model,Removed_term,R2,Adjusted_R2,RMSE,' ...
        'Mean_residual,Residual_std,Durbin_Watson,Max_abs_ACF,' ...
        'ACF_95_limit,Ljung_Box_Q,Ljung_Box_p,Jarque_Bera,' ...
        'Jarque_Bera_p,Max_abs_regressor_residual_corr,AIC,BIC\n']);
    for outputIndex = 1:numel(responseNames)
        modelNames = {'Initial_OLS','Final_RLS','Reduced'};
        fits = {OLS(outputIndex),RLS(outputIndex),Reduced(outputIndex)};
        removed = {'','',parameterNames{leastInfluentialIndex(outputIndex)}};
        for modelIndex = 1:3
            m = fits{modelIndex}.validation;
            fprintf(fileID,['%s,%s,%s,%.12g,%.12g,%.12g,%.12g,' ...
                '%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,' ...
                '%.12g,%.12g,%.12g\n'], ...
                responseNames{outputIndex},modelNames{modelIndex}, ...
                removed{modelIndex},m.R2,m.adjustedR2,m.RMSE, ...
                m.meanResidual,m.standardDeviationResidual,m.DurbinWatson, ...
                m.maxAbsoluteAutocorrelation,m.autocorrelation95Limit, ...
                m.LjungBoxQ,m.LjungBoxPValue,m.JarqueBera, ...
                m.JarqueBeraPValue, ...
                m.maxAbsoluteRegressorResidualCorrelation,m.AIC,m.BIC);
        end
    end
    fclose(fileID);
end

function writeDominanceCSV(filePath,responseNames,parameterNames, ...
        contribution,fraction)
    fileID = fopen(filePath,'w');
    if fileID < 0; error('Could not create %s.',filePath); end
    fprintf(fileID,'Response,Parameter,Contribution,Dynamic_fraction\n');
    for outputIndex = 1:numel(responseNames)
        for parameterIndex = 1:numel(parameterNames)
            fprintf(fileID,'%s,%s,%.12g,%.12g\n', ...
                responseNames{outputIndex},parameterNames{parameterIndex}, ...
                contribution(parameterIndex,outputIndex), ...
                fraction(parameterIndex,outputIndex));
        end
    end
    fclose(fileID);
end

function writeIdentifiabilityCSV(filePath,parameterNames,responseNames,RLS, ...
        conditionNumber,VIF,correlationMatrix)
    fileID = fopen(filePath,'w');
    if fileID < 0; error('Could not create %s.',filePath); end
    fprintf(fileID,'Metric,Term_1,Term_2,Value\n');
    fprintf(fileID,'Normalized_condition_number,,,%.12g\n',conditionNumber);
    for index = 2:numel(parameterNames)
        fprintf(fileID,'VIF,%s,,%.12g\n',parameterNames{index},VIF(index-1));
    end
    for row = 1:size(correlationMatrix,1)
        for column = row+1:size(correlationMatrix,2)
            fprintf(fileID,'Regressor_correlation,%s,%s,%.12g\n', ...
                parameterNames{row+1},parameterNames{column+1}, ...
                correlationMatrix(row,column));
        end
    end
    for outputIndex = 1:numel(responseNames)
        parameterCorrelation = RLS(outputIndex).parameterCorrelation;
        for row = 1:size(parameterCorrelation,1)
            for column = row+1:size(parameterCorrelation,2)
                fprintf(fileID,'%s_parameter_correlation,%s,%s,%.12g\n', ...
                    responseNames{outputIndex},parameterNames{row}, ...
                    parameterNames{column},parameterCorrelation(row,column));
            end
        end
        certainty = abs(RLS(outputIndex).theta(2:end)) ...
            ./max(RLS(outputIndex).hacStandardError(2:end),eps);
        [mostCertainValue,mostCertainIndex] = max(certainty);
        [leastCertainValue,leastCertainIndex] = min(certainty);
        fprintf(fileID,'%s_most_certain,%s,,%.12g\n', ...
            responseNames{outputIndex}, ...
            parameterNames{mostCertainIndex+1},mostCertainValue);
        fprintf(fileID,'%s_least_certain,%s,,%.12g\n', ...
            responseNames{outputIndex}, ...
            parameterNames{leastCertainIndex+1},leastCertainValue);
    end
    fclose(fileID);
end

function writeDataUseCSV(filePath,rawFile,ekfFile,t,initialMask, ...
        updateMask,validationMask)
    fileID = fopen(filePath,'w');
    if fileID < 0; error('Could not create %s.',filePath); end
    fprintf(fileID,'Item,Source_or_block,Purpose,First_time_s,Last_time_s,Samples\n');
    fprintf(fileID,'File,%s,Commanded controls and timing,,,\n',rawFile);
    fprintf(fileID,['File,%s,EKF states corrected IMU and estimated ' ...
        'biases,,,\n'],ekfFile);
    masks = {initialMask,updateMask,validationMask};
    names = {'Initial_OLS','RLS_update','Independent_validation'};
    purposes = {'Initial parameter estimate', ...
        'Chronological recursive update', ...
        'Unseen model validation'};
    for index = 1:3
        mask = masks{index};
        fprintf(fileID,'Block,%s,%s,%.6f,%.6f,%d\n', ...
            names{index},purposes{index},min(t(mask)),max(t(mask)),nnz(mask));
    end
    fclose(fileID);
end
