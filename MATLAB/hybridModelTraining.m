%% General preparations 
clear 
close all
clc

%% Generation of database for training
generateDatabase;
V = 0.001;

% Split data into training, validation, and testing sets
XTrain = X(trainInd,:);
YTrain = Y(trainInd,:);
XVal = X(valInd,:);
YVal = Y(valInd,:);
XTest = X(testInd:end,:);
YTest = Y(testInd:end,:);

% Downsampling data for training
[XTrainDownsampled, YTrainDownsampled] = downsampleData(XTrain, YTrain, 5);

% Normalization of data
XTrainDownsampled = cell2mat(XTrainDownsampled);
YTrainDownsampled = cell2mat(YTrainDownsampled);

[XTrainDownsampledNormalized, YTrainDownsampledNormalized] = normalizeData(XTrainDownsampled, YTrainDownsampled);

%% Creation of ANN
% Three hidden layers with 8 neurons each network
hiddenLayerSize = [8];
netAceticAcid = feedforwardnet(hiddenLayerSize, 'trainbfg');
netMethanol = feedforwardnet(hiddenLayerSize, 'trainbfg');
netMethylAcetate = feedforwardnet(hiddenLayerSize, 'trainbfg');
netWater = feedforwardnet(hiddenLayerSize, 'trainbfg');

% Set activation function as sigmoid for each layer
for i = 1:width(hiddenLayerSize)
    netAceticAcid.layers{i}.transferFcn = 'logsig';
    netMethanol.layers{i}.transferFcn = 'logsig';
    netMethylAcetate.layers{i}.transferFcn = 'logsig';
    netWater.layers{i}.transferFcn = 'logsig';
end

% Training
netAceticAcid = train(netAceticAcid, XTrainDownsampledNormalized', YTrainDownsampledNormalized(:, 1)');
netMethanol = train(netMethanol, XTrainDownsampledNormalized', YTrainDownsampledNormalized(:, 2)');
netMethylAcetate = train(netMethylAcetate, XTrainDownsampledNormalized', YTrainDownsampledNormalized(:, 3)');
netWater = train(netWater, XTrainDownsampledNormalized', YTrainDownsampledNormalized(:, 4)');

% Save the trained networks
save('trained_ANNs.mat', 'netAceticAcid', 'netMethanol', 'netMethylAcetate', 'netWater');

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
