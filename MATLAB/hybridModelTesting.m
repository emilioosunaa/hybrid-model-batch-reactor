%% General preparations 
clear 
close all
clc

%% Trainning of the hybrid model
if ~exist('netAceticAcid', 'var') || ~exist('netMethanol', 'var') || ...
   ~exist('netMethylAcetate', 'var') || ~exist('netWater', 'var')
    hybridModelTraining;

else
    fprintf('Neural networks already exist. Skipping training.\n');
end

%% Generation of database for testing
XTest = cell2mat(XTest);
YTest = cell2mat(YTest);

% Normalization of data
[XTestNormalized, YTestNormalized] = normalizeData(XTest, YTest);

%% Testing of the hybrid model
% Making predictions
YPredAceticAcidNormalized = netAceticAcid(XTestNormalized');
YPredMethanolNormalized = netMethanol(XTestNormalized');
YPredMethylAcetateNormalized = netMethylAcetate(XTestNormalized');
YPredWaterNormalized = netWater(XTestNormalized');

YPredNormalized = [YPredAceticAcidNormalized; YPredMethanolNormalized; YPredMethylAcetateNormalized; YPredWaterNormalized]';

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
minValueY = [0 0 0 0];
maxValueY = [20/V 20/V 20/V 20/V];
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