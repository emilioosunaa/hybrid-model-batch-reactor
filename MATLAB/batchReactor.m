%% ODE system: batch reactor
% Esterification reaction of methanol by acetic acid.
function dy = batchReactor(t, y, T, QCat)
    % Kinetic parameters
    k0e = 4.21;         % [m^3 mol^-1 s^-1 mL_cat^-1]
    k0h = 0.322;        % [m^3 mol^-1 s^-1 mL_cat^-1]
    Eae = 53804.0;      % [J mol^-1]
    Eah = 52584.0;      % [J mol^-1]
    
    % Constants
    R = 8.314;          % [J K^-1 mol]

    % Unpacking of temperature [K], amount of catalyst [mL],
    % volume [m^3] and concentrations [mol m^-3]
    CAA = y(1);
    CM = y(2);
    CMA = y(3);
    CW = y(4);

    % Reaction rate [mol s^-1 m^-3]
    r = k0e*QCat*exp(-Eae/R/T)*CAA*CM - k0h*QCat*exp(-Eah/R/T)*CMA*CW;
    
    % Mass balance equations
    dy = zeros(4, 1);
    dy(1) = -r;       % dCAA_dt
    dy(2) = -r;       % dCM_dt
    dy(3) = r;        % dCMA_dt
    dy(4) = r;        % dCW_dt
end