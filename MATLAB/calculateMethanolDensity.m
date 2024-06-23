function rhoM = calculateMethanolDensity(T)
    % Constants and formula from Perrys Chemical Engineers Handbook 8th ed.
    C1 = 2.3267;
    C2 = 0.27073;
    C3 = 512.5;
    C4 = 0.24713;
  
    rhoM = C1/(C2^(1 + (1-T/C3)^C4)) * 1000;   % [mol/m^3]
end