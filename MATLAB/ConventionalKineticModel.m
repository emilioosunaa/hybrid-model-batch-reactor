%% General preparations 
clear 
close all
clc

%%  Conventional Kinetic Model
% Define the initial conditions
% T[°C] Qcat[mL] CAA[mol L^-1] CW[mol L^-1] feeding time[s] FM[mol s^-1] Number of moles of M fed [mol]
initialConditions = [
    66.55 9.14 16.06 3.57 200 0.0539 10.79;
    80.00 7.50 14.86 7.75 200 0.0288 5.77;
    60.00 8.00 13.15 13.15 200 0.0305 6.11;
    60.00 8.00 13.15 13.15 600 0.0101 6.11
];

% Operation conditions
CAA0 = initialConditions(:, 3) * 1000;      % [mol m^-3]
CM0 = 0.00;                                 % [mol m^-3]
CMA0 = 0.00;                                % [mol m^-3]
CW0 = initialConditions(:, 4) * 1000;       % [mol m^-3]
QCat = initialConditions(:, 2);             % [mL]
T = initialConditions(:, 1) + 273.15;       % [K]
FM = initialConditions(:, 6);               % [mol s^-1]
V0 = 0.001;                                 % [m^-3]

% Time settings
dt = 1.0;                                   % [s]
feedingTime = initialConditions(1, 5);      % [s]
totalTime = 2000;                           % [s]
tspan1 = [0, feedingTime];
tspan2 = [feedingTime, totalTime];

y0 = [CAA0(1)*V0, CM0*V0, CMA0*V0, CW0(1)*V0, V0];

% Solve the differential equations during feed
[t1, y1] = ode45(@fedbatchReactor, tspan1, y0, [], T(1), QCat(1), FM(1));

% Get final conditions from first phase as initial conditions for second phase
y0AfterFeed = [y1(end, 1)/y1(end, 5), y1(end, 2)/y1(end, 5), y1(end, 3)/y1(end, 5), y1(end, 4)/y1(end, 5)];

% Solve the differential equations after feed is stopped
[t2, y2] = ode45(@batchReactor, tspan2, y0AfterFeed, [], T(1), QCat(1));

% Combine results
t = [t1; t2];
y = [y1(:, 1:4)./y1(:, 5); y2];
tDownsampled = t(1:5:end, :);
yDownsampled = y(1:5:end, :);
% Plot the results
figure;
hold all;
plot(tDownsampled, yDownsampled(:,1)/1000, '-+', 'MarkerSize', 4);
plot(tDownsampled, yDownsampled(:,2)/1000, '-x', 'MarkerSize', 4);
plot(tDownsampled, yDownsampled(:,3)/1000, '-*', 'MarkerSize', 4);
plot(tDownsampled, yDownsampled(:,4)/1000, '-o', 'MarkerSize', 4);
legend('C AA', 'C M', 'C MA', 'C W');
xlabel('Time (s)');
ylabel('Concentration [mol L^-1]');
title('Fed-Batch Reactor Simulation');
grid on;
hold off;
