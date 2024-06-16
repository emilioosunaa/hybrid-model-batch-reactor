%% General preparations 
clear 
close all
clc


%% Generation of database for training, testing and validation
% Extraction of initial conditions for the runs
initialConditions = readmatrix('initial_conditions.csv');

V = 0.001;   % [m^3]
[X, Y] = generateDatabase(initialConditions, V);

% Split data into training, validation, and testing sets
trainInd = 1:16;
valInd = 17:21;
testInd = 22;

XTrain = X(trainInd,:);
YTrain = Y(trainInd,:);
XVal = X(valInd,:);
YVal = Y(valInd,:);
XTest = X(testInd:end,:);
YTest = Y(testInd:end,:);

% Downsampling data for training
downsampleFactor = 5;
[XTrainDownsampled, YTrainDownsampled] = downsampleData(XTrain, YTrain, downsampleFactor);

% Normalization of data
XTrainDownsampled = cell2mat(XTrainDownsampled);
YTrainDownsampled = cell2mat(YTrainDownsampled);
XTest = cell2mat(XTest);
YTest = cell2mat(YTest);
XVal = cell2mat(XVal);
YVal = cell2mat(YVal);

[XTrainDownsampledNormalized, YTrainDownsampledNormalized] = normalizeData(XTrainDownsampled, YTrainDownsampled, V);
[XTestNormalized, YTestNormalized] = normalizeData(XTest, YTest, V);
[XValNormalized, YValNormalized] = normalizeData(XVal, YVal, V);

%% Creation of ANN
% Create a dictionary (containers.Map) of neural networks
speciesNames = {'Acetic Acid', 'Methanol', 'Methyl Acetate', 'Water'};
networks = containers.Map();

% Initialize lists for the predictions
YTestPredNormalized = zeros(size(YTestNormalized));
YValPredNormalized = zeros(size(YValNormalized));

% Three hidden layers with 8 neurons each network, activation function as sigmoid for each layer
hiddenLayerSize = [8];
for i = 1:length(speciesNames)
    name = speciesNames{i};
    net = feedforwardnet(hiddenLayerSize, 'trainbfg');

    for j = 1:width(hiddenLayerSize)
        net.layers{j}.transferFcn = 'logsig';
    end

    % Setup Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 100/100;
    net.divideParam.valRatio = 0/100;
    net.divideParam.testRatio = 0/100;

    % Train the network
    [net, tr] = train(net, XTrainDownsampledNormalized', YTrainDownsampledNormalized(:, i)');

    % Validate the network
    YValPredNormalized(:, i) = net(XValNormalized');
    %perfVal = perform(net, YValNormalized(:, i)'', YValPredNormalized(:, i)');
    %fprintf('Validation Performance for %s: %f\n', name, perfVal);
    
    % Test the network
    YTestPredNormalized(:, i) = net(XTestNormalized');
    %perfTest = perform(net, YTestNormalized(:, i)'', YTestPredNormalized(:, i)');
    %fprintf('Test Performance for %s: %f\n', name, perfTest);

    % Save the network in the dictionary
    networks(name) = net;
end

%% Calculate evaluation metrics
% Initialize variables
nSamples = 1000;       % Number of samples
mseTest = 0;
mseVal = 0;

% Calculate mean square error for the testing part
for j = 1:4
    for i = 1:nSamples
        mseTest = mseTest + (YTestPredNormalized(i, j) - YTestNormalized(i, j))^2;
        mseVal = mseVal + (YValPredNormalized(i, j) - YValNormalized(i, j))^2;
    end
end

mseTest = mseTest / (4 * nSamples);
mseVal = mseVal / (4 * nSamples);

fprintf('Mean Squared Error from testing: %.6e\n', mseTest);
fprintf('Mean Squared Error from validation: %.6e\n', mseVal);

% Plot predictions
plotPredictions(YTestPredNormalized, XTest, YTest, 16, V, 50);

% Save the trained networks
save('trained_networks.mat', 'networks');



%% Helper functions
function [XDownsampled, YDownsampled] = downsampleData(X, Y, downsampleFactor)
    columns = width(X);
    rows = height(X);
    XDownsampled = cell(1, columns);
    YDownsampled = cell(1, 4);
    
    for i = 1:rows
         % Downsampling Y
        for col = 4:columns
            currentVector = Y{i, col - 3};
    
            % Downsample the vector by selecting every {downsampleFactor}th element
            downsampledVector = currentVector(1:downsampleFactor:end);
    
            % Store the downsampled vector back in the corresponding position
            YDownsampled{i, col - 3} = downsampledVector;
        end
        
        % Downsampling X
        for col = 1:columns
            % Extract the current vector
            currentVector = X{i, col};
    
            % Downsample the vector by selecting every 5th element
            downsampledVector = currentVector(1:downsampleFactor:end);
    
            % Store the downsampled vector back in the corresponding position
            XDownsampled{i, col} = downsampledVector;
        end

    end

end
