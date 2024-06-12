package HybridModel

  model BatchReactor
    // Parameters
    parameter Real T = 298.15 "Temperature [K]";
    parameter Real QCat = 1 "Amount of catalyst [mL]";
    parameter Real R = 8.314 "Universal gas constant [J mol^-1 K^-1]";
    parameter Real k0e = 1e3 "Pre-exponential factor for esterification [m^3 mol^-1 s^-1]";
    parameter Real k0h = 1e3 "Pre-exponential factor for hydrolysis [m^3 mol^-1 s^-1]";
    parameter Real Eae = 50000 "Activation energy for esterification [J mol^-1]";
    parameter Real Eah = 60000 "Activation energy for hydrolysis [J mol^-1]";
  
    // State variables
    Real CAA(start=1) "Concentration of acetic acid [mol m^-3]";
    Real CM(start=1) "Concentration of methanol [mol m^-3]";
    Real CMA(start=0) "Concentration of methyl acetate [mol m^-3]";
    Real CW(start=0) "Concentration of water [mol m^-3]";
  
    // Reaction rate [mol s^-1 m^-3]
    Real r;
  
  equation
    // Define the reaction rate
    r = k0e * QCat * exp(-Eae / (R * T)) * CAA * CM - k0h * QCat * exp(-Eah / (R * T)) * CMA * CW;
  
    // Mass balance equations
    der(CAA) = -r;
    der(CM) = -r;
    der(CMA) = r;
    der(CW) = r;
  
  end BatchReactor;


end HybridModel;
