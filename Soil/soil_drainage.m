function [f1] = soil_drainage(wa_unsat, theta_sat, ks, Dmin, Dmax)
    % Ducharne % Polcher, 1998
    dd = 1.5;
    %% layer1
    % Dmin = 0.048; % mm day-1
    % Dmax = 4.8; % mm day-1
    % 
    %% layer2
    % Dmin = 0.012; % 0.0005*24, mm day-1
    % Dmax = 1.2; % 0.05*24,   mm day-1
    thx1 = wa_unsat / theta_sat;
    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end
    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);
end
