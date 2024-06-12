%% General preparations 
%clear 
close all
clc

%% Trainning of the hybrid model
if ~exist('netAceticAcid', 'var') || ~exist('netMethanol', 'var') || ...
   ~exist('netMethylAcetate', 'var') || ~exist('netWater', 'var')
    hybridmodel_training;

else
    fprintf('Neural networks already exist. Skipping training.\n');
end

%% Generation of database for testing
rows = height(initialConditionsTest);

% Operation conditions
V = 0.001;                                          % [m^3]
dt = 2000;                                          % [s] 
CAA0 = initialConditionsTest(:, 1) ./ V;           % [mol]
CM0 = initialConditionsTest(:, 2) ./ V;            % [mol]
CMA0 = 0.00;                                        % [mol]
CW0 = initialConditionsTest(:, 3) ./ V;            % [mol]
QCat = initialConditionsTest(:, 4);                % [mL]
T = initialConditionsTest(:, 5) + 273.15;          % [K]

% Initialization of the data base
columns = 7;
XTest = cell(1, columns);
YTest = cell(1, 4);
tspan = linspace(0, dt, 1000);

for i = 1:rows
    % Initial conditions (C_AA, C_M, C_MA, C_W)
    y0 = [CAA0(i), CM0(i), CMA0, CW0(i)];

    % ODE solution
    [t, y] = ode45(@batchReactor, tspan, y0, [], T(i), QCat(i));
    
    % Saving the data
    XTest{i, 1} = t;
    XTest{i, 2} = zeros(size(t)) + T(i);
    XTest{i, 3} = zeros(size(t)) + QCat(i);
    XTest{i, 4} = y(:, 1);
    XTest{i, 5} = y(:, 2);
    XTest{i, 6} = y(:, 3);
    XTest{i, 7} = y(:, 4);

    % Rolling XTest so t+1 data for YTest is obtained
    for col = 4:columns
        % Extract the current data for the specific experiment and column
        currentData = XTest{i, col};

        % Shift data, assuming there is enough data
        YTest{i, col - 3} = [currentData(2:end); currentData(end)];
    end
end

XTest = cell2mat(XTest);
YTest = cell2mat(YTest);

% Normalization of data
minValueX = [0 273 0 0 0 0 0];
maxValueX = [1 373 15 20/V 20/V 20/V 20/V];
XTestNormalized = (XTest - minValueX) ./ (maxValueX - minValueX);

%% Testing of the hybrid model
% Making predictions
YPredAceticAcidNormalized = netAceticAcid(XTestNormalized');
YPredMethanolNormalized = netMethanol(XTestNormalized');
YPredMethylAcetateNormalized = netMethylAcetate(XTestNormalized');
YPredWaterNormalized = netWater(XTestNormalized');

YPredNormalized = [YPredAceticAcidNormalized; YPredMethanolNormalized; YPredMethylAcetateNormalized; YPredWaterNormalized]';

minValueY = [0 0 0 0];
maxValueY = [20/V 20/V 20/V 20/V];
YTestNormalized = (YTest - minValueY) ./ (maxValueY - minValueY);

%% Calculate evaluation metrics
% Initialize MSE
nSamples = 1000;       % Number of samples
mse = 0;

% Calculate MSE
for j = 1:4
    for i = 1:nSamples
        mse = mse + (YPredNormalized(i, j) - YTestNormalized(i, j))^2;
    end
end

mse = mse / (4 * nSamples);

fprintf('Mean Squared Error: %.6e\n', mse);

%% Plotting predictions
downsampleFactor = 50;
YPred = YPredNormalized .* (maxValueY - minValueY) + minValueY;
YTestDownsampled = YTest(15001:downsampleFactor:16000, :)*0.001;
YPredDownsampled = YPred(15001:downsampleFactor:16000, :)*0.001;
timeDownsampled = XTest(15001:downsampleFactor:16000, 1);

figure;
hold all

for i = 1:4
    plot(timeDownsampled, YTestDownsampled(:, i), 'o', 'MarkerSize', 4);
    plot(timeDownsampled, YPredDownsampled(:, i), 'LineWidth', 2);
end

title ('Esterification reaction of methanol by acetic acid in a batch reactor');
legend('C AA', 'C M', 'C MA', 'C W');
xlabel('Time [s]');
ylabel('Concentration [mol/L]');
grid on;
hold off;