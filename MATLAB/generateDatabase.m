% Extraction of initial conditions for the runs
initialConditions = readmatrix('initial_conditions.csv');

%% Generation of database for training
trainInd = 1:16;
valInd = 17:21;
testInd = 22;
rows = height(initialConditions);

% Operation conditions
V = 0.001;                                     % [m^3]
dt = 2000;                                     % [s] 
CAA0 = initialConditions(:, 1) ./ V;           % [mol]
CM0 = initialConditions(:, 2) ./ V;            % [mol]
CMA0 = 0.00;                                   % [mol]
CW0 = initialConditions(:, 3) ./ V;            % [mol]
QCat = initialConditions(:, 4);                % [mL]
T = initialConditions(:, 5) + 273.15;          % [K]

% Initialization of the data base
columns = 7;
X = cell(1, columns);
Y = cell(1, 4);
tspan = linspace(0, dt, 1000);

for i = 1:rows
    % Initial conditions (C_AA, C_M, C_MA, C_W)
    y0 = [CAA0(i), CM0(i), CMA0, CW0(i)];

    % ODE solution
    [t, y] = ode45(@batchReactor, tspan, y0, [], T(i), QCat(i));
    
    % Saving the data
    X{i, 1} = t;
    X{i, 2} = zeros(size(t)) + T(i);
    X{i, 3} = zeros(size(t)) + QCat(i);
    X{i, 4} = y(:, 1);
    X{i, 5} = y(:, 2);
    X{i, 6} = y(:, 3);
    X{i, 7} = y(:, 4);

    % Rolling X so t+1 data for Y is obtained
    for col = 4:columns
        % Extract the current data for the specific experiment and column
        currentData = X{i, col};

        % Shift data, assuming there is enough data
        Y{i, col - 3} = [currentData(2:end); currentData(end)];
    end

end