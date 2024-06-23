function dy = fedbatchReactor(t, y, T, QCat, FM)    
    % Kinetic parameters
    k0e = 4.21;         % [m^3 mol^-1 s^-1 mL_cat^-1]
    k0h = 0.322;        % [m^3 mol^-1 s^-1 mL_cat^-1]
    Eae = 53804.0;      % [J mol^-1]
    Eah = 52584.0;      % [J mol^-1]
    
    % Constants
    R = 8.314;                                      % [J K^-1 mol]
    MM = 32.042/1000;                               % [kg mol^-1]
    rhoM = calculateMethanolDensity(T) * MM;        % [kg/m^3]
    
    % Extract variables
    NAA = y(1);         % mol
    NM = y(2);          % mol
    NMA = y(3);         % mol
    NW = y(4);          % mol
    V = y(5);           % m^3
    
    % Calculate Concentrations
    CAA = NAA / V;
    CM = NM / V;
    CMA = NMA / V;
    CW = NW / V;

    % Reaction rate [mol s^-1 m^-3]
    r = k0e*QCat*exp(-Eae/R/T)*CAA*CM - k0h*QCat*exp(-Eah/R/T)*CMA*CW;
    
    % Differential equations
    dy = zeros(5, 1);
    dy(1) = -V * r;         % dNAA_dt
    dy(2) = FM - V * r;     % dNM_dt
    dy(3) = V * r;          % dNMA_dt
    dy(4) = V * r;          % dNW_dt
    Q_t = FM * MM / rhoM;   % Q_t
    dy(5) = Q_t;            % dV_dt
end