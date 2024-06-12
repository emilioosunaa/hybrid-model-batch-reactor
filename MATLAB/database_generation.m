%% General preparations 
clear 
close all
clc

% Extraction of initial conditions for the runs
initialConditionsFileName = 'initial_conditions.csv';
initialConditions = readmatrix(initialConditionsFileName);

% Kinetic parameters
k0e = 4.21;         % [m^3 mol^-1 s^-1 mL_cat^-1]
k0h = 0.322;        % [m^3 mol^-1 s^-1 mL_cat^-1]
Eae = 53804.0;      % [J mol^-1]
Eah = 52584.0;      % [J mol^-1]

% Constants
R = 8.314;          % [J K^-1 mol]

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
tspan = linspace(0, dt, 1000);

for i = 1:rows
    % Initial conditions (C_AA, C_M, C_MA, C_W)
    y0 = [CAA0(i), CM0(i), CMA0, CW0(i)];

    % ODE solution
    [t, y] = ode45(@batchReactor, tspan, y0, [], T(i), QCat(i), R, k0e, k0h, Eae, Eah);
    
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
        downsampledVector = currentVector(1:5:end);

        % Store the downsampled vector back in the corresponding position
        XTrainDownsampled{i, col} = downsampledVector;
    end

end

% Plotting training data
figure;
hold all
plot(XTrainDownsampled{1,1}, XTrainDownsampled{1,4}*0.001, 'LineWidth', 2);
plot(XTrainDownsampled{1,1}, XTrainDownsampled{1,5}*0.001, 'LineWidth', 2);
plot(XTrainDownsampled{1,1}, XTrainDownsampled{1,6}*0.001, 'LineWidth', 2);
plot(XTrainDownsampled{1,1}, XTrainDownsampled{1,7}*0.001, 'LineWidth', 2);
title ('Esterification reaction of methanol by acetic acid in a batch reactor');
legend('AA','M','MA','W');
xlabel('Time [s]');   ylabel('Concentration [mol/L]');

XTrainDownsampled = cell2mat(XTrainDownsampled);
YTrainDownsampled = cell2mat(YTrainDownsampled);

% Normalization of data
minValueX = [0 273 0 0 0 0 0];
maxValueX = [1 373 15 20/V 20/V 20/V 20/V];
XTrainDownsampledNormalized = (XTrainDownsampled - minValueX) ./ (maxValueX - minValueX);

minValueY = [0 0 0 0];
maxValueY = [20/V 20/V 20/V 20/V];
YTrainDownsampledNormalized = (YTrainDownsampled - minValueY) ./ (maxValueY - minValueY);
