%% General preparations 
clear 
close all
clc

% Extraction of initial conditions for the runs
initialConditionsFileName = 'initial_conditions.csv';
initialConditions = readmatrix(initialConditionsFileName);
initialConditionsTrain = initialConditions(1:16,:);
initialConditionsValidation = initialConditions(17:21,:);
initialConditionsTest = initialConditions(22:end,:);

%% Generation of database for training
rows = height(initialConditionsTrain);

% Operation conditions
V = 0.001;                                          % [m^3]
dt = 2000;                                          % [s] 
CAA0 = initialConditionsTrain(:, 1) ./ V;           % [mol]
CM0 = initialConditionsTrain(:, 2) ./ V;            % [mol]
CMA0 = 0.00;                                        % [mol]
CW0 = initialConditionsTrain(:, 3) ./ V;            % [mol]
QCat = initialConditionsTrain(:, 4);                % [mL]
T = initialConditionsTrain(:, 5) + 273.15;          % [K]

% Initialization of the database
columns = 7;
downsampleFactor = 5;
XTrain = cell(1, columns);
YTrain = cell(1, 4);
XTrainDownsampled = cell(size(XTrain));
YTrainDownsampled = cell(size(YTrain));
tspan = linspace(0, dt, 1000);

for i = 1:rows
    % Initial conditions (C_AA, C_M, C_MA, C_W)
    y0 = [CAA0(i), CM0(i), CMA0, CW0(i)];

    % ODE solution
    [t, y] = ode45(@batchReactor, tspan, y0, [], T(i), QCat(i));
    
    % Saving the data
    XTrain{i, 1} = t;
    XTrain{i, 2} = zeros(size(t)) + T(i);
    XTrain{i, 3} = zeros(size(t)) + QCat(i);
    XTrain{i, 4} = y(:, 1);
    XTrain{i, 5} = y(:, 2);
    XTrain{i, 6} = y(:, 3);
    XTrain{i, 7} = y(:, 4);
    
    % Rolling XTrain so t+1 data for YTrain is obtained
    for col = 4:columns
        % Extract the current data for the specific experiment and column
        currentData = XTrain{i, col};

        % Shift data, assuming there is enough data
        YTrain{i, col - 3} = [currentData(2:end); currentData(end)];
        
        % Downsampling
        currentVector = YTrain{i, col - 3};

        % Downsample the vector by selecting every {downsampleFactor}th element
        downsampledVector = currentVector(1:downsampleFactor:end);

        % Store the downsampled vector back in the corresponding position
        YTrainDownsampled{i, col - 3} = downsampledVector;
    end
    
    % Downsampling
    for col = 1:columns
        % Extract the current vector
        currentVector = XTrain{i, col};

        % Downsample the vector by selecting every 5th element
        downsampledVector = currentVector(1:downsampleFactor:end);

        % Store the downsampled vector back in the corresponding position
        XTrainDownsampled{i, col} = downsampledVector;
    end

end

XTrainDownsampled = cell2mat(XTrainDownsampled);
YTrainDownsampled = cell2mat(YTrainDownsampled);

% Normalization of data
minValueX = [0 273 0 0 0 0 0];
maxValueX = [1 373 15 20/V 20/V 20/V 20/V];
XTrainDownsampledNormalized = (XTrainDownsampled - minValueX) ./ (maxValueX - minValueX);

minValueY = [0 0 0 0];
maxValueY = [20/V 20/V 20/V 20/V];
YTrainDownsampledNormalized = (YTrainDownsampled - minValueY) ./ (maxValueY - minValueY);

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
