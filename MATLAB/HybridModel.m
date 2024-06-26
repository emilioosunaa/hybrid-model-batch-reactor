%% General preparations 
clear 
close all
clc

%%  Hybrid Model
% Load the pre-trained ANN models
load('trained_networks.mat'); % Load trained ANN, that are stored in a dictionary (containers.Map) of neural networks 

%concentrationsTest = calculateConcentrationChanges(networks, XTest);

% Define the initial conditions
% T[°C] Qcat[mL] CAA[mol L^-1] CW[mol L^-1] feeding time[s] FM[mol s^-1] Number of moles of M fed [mol]
initialConditions = [
    66.55 9.14 16.06 3.57 200 0.0539 10.79;
    80.00 7.50 14.86 7.75 200 0.0288 5.77;
    60.00 8.00 13.15 13.15 200 0.0305 6.11;
    60.00 8.00 13.15 13.15 600 0.0101 6.11
];

% Constansts
MM = 32.042/1000;                           % [kg mol^-1], methanol molecular weight

% Operation conditions
CAA0 = initialConditions(:, 3) * 1000;      % [mol m^-3]
CM0 = 0.00;                                 % [mol m^-3]
CMA0 = 0.00;                                % [mol m^-3]
CW0 = initialConditions(:, 4) * 1000;       % [mol m^-3]
QCat0 = initialConditions(:, 2);            % [mL]
T0 = initialConditions(:, 1) + 273.15;      % [K]
FM0 = initialConditions(:, 6);              % [mol s^-1]
V0 = 0.001;                                 % [m^-3]

% Time settings
dt = 0.5;                                   % [s]
feedingTime = initialConditions(1, 5);      % [s] CHANGE THIS
totalTime = 2000;                           % [s]
numStepsFeed = feedingTime / dt;
numStepsTotal = totalTime / dt;

% Preallocate arrays for results
time = (0:dt:totalTime)';
NAA = zeros(numStepsTotal + 1, 1);
NM = zeros(numStepsTotal + 1, 1);
NMA = zeros(numStepsTotal + 1, 1);
NW = zeros(numStepsTotal + 1, 1);
V = zeros(numStepsTotal + 1, 1);

% Initial values
T = T0(1);
QCat = QCat0(1);
NAA(1) = CAA0(1) * V0;
NM(1) = CM0 * V0;
NMA(1) = CMA0 * V0;
NW(1) = CW0(1) * V0;
V(1) = V0;
FM = FM0(1);
rhoM = calculateMethanolDensity(T) * MM; % [kg/m^3]

% Function to compute concentration changes using ANN models
function concentrations = calculateConcentrationChanges(networks, inputs)
    % Extracting networks
    netAA = networks('Acetic Acid');
    netM = networks('Methanol');
    netMA = networks('Methyl Acetate');
    netW = networks('Water');
    
    % Normalizing inputs
    minValueInput = [0 273 0 0 0 0 0];
    maxValueInput = [1 373 15 20/0.001 20/0.001 20/0.001 20/0.001];
    inputsNormalized = (inputs - minValueInput) ./ (maxValueInput - minValueInput);

    % Predicting concentrations for t+1, results from ANN are normalized.
    CAA = netAA(inputsNormalized');
    CM = netM(inputsNormalized');
    CMA = netMA(inputsNormalized');
    CW = netW(inputsNormalized');

    concentrationsNormalized = [CAA' CM' CMA' CW'];

    % Denormalizing outputs
    minValueOutput = [0 0 0 0];
    maxValueOutput = [20/0.001 20/0.001 20/0.001 20/0.001];
    concentrations = concentrationsNormalized .* (maxValueOutput - minValueOutput) + minValueOutput;
end

% Simulation loop
for t = 1:numStepsTotal
    % Current concentrations
    CAA_t = NAA(t) / V(t);
    CM_t = NM(t) / V(t);
    CMA_t = NMA(t) / V(t);
    CW_t = NW(t) / V(t);
    
    % Inputs for ANN models
    inputs = [t*dt T QCat CAA_t CM_t CMA_t CW_t];
    
    % Step 1: Compute concentration changes using ANN models
    C_i0_t1 = calculateConcentrationChanges(networks, inputs);

    % Step 2: Add the feeding term (only during feeding time)
    if t <= numStepsFeed
        NM(t + 1) = C_i0_t1(2) * V(t) + FM * dt;
    else
        NM(t + 1) = C_i0_t1(2) * V(t);
    end

    NAA(t + 1) = C_i0_t1(1) * V(t);
    NMA(t + 1) = C_i0_t1(3) * V(t);
    NW(t + 1) = C_i0_t1(4) * V(t);

    % Step 3: Correction of the volume
    if t <= numStepsFeed
        V(t + 1) = V(t) + FM * MM / rhoM * dt;
    else
        V(t + 1) = V(t);
    end
    
end

% Compute final concentrations
CAA = NAA ./ V;
CM = NM ./ V;
CMA = NMA ./ V;
CW = NW ./ V;

% Plot the results
figure;
hold all;
plot(time, CAA/1000, '-', 'MarkerSize', 4);
plot(time, CM/1000, '-', 'MarkerSize', 4);
plot(time, CMA/1000, '-', 'MarkerSize', 4);
plot(time, CW/1000, '-', 'MarkerSize', 4);
legend('C AA', 'C M', 'C MA', 'C W');
xlabel('Time (s)');
ylabel('Concentration [mol L^-1]');
title('Hybrid Model Fed-Batch Reactor Simulation');
grid on;
hold off;