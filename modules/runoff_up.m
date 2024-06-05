function [srf, IWS, Vmax] = runoff_up(Pnet, zgw, zm, wa, soilpar)
%%INPUT:
% zm       : Soil depth of different layers (mm)
% wa       : The antecedent soil water content (mm) expressed
%            as a function of the WHC in that layer
% soilpar  : Soil parameters according to Soil type
% Pnet     : Net precipitation = P-I+Snowmelt
% zgw      : groundwater table depth
%
%%OUTPUT:
% srf      : Surface Runoff, mm
% IWS      : Water enter into soil surface, mm
% 
%%Reference:
% 1. Choudhury BJ, Digirolamo NE, 1998
%
% 2. SCS (1985). National Engineering Handbook. Section 4: Hydrology.
% Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.
%
% 3. SCS (1986). Urban Hydrology for Small Watersheds, Technical Release No. 55.
% Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.

theta_sat = soilpar(3); % saturated wa

if zgw <= 0
    % exceeded groundwater on the soil surface
    Vmax = 0;
    srf = 0 - zgw * theta_sat;
    Pnet = 0;
else
    if zgw > 0 && zgw <= zm(1)
        d1 = zgw; % the thickness of unsaturated soil, (mm)
        
        % the unsaturated soil water, (mm)
        wa1_unsat = (wa(1) * zm(1) - theta_sat * (zm(1) - d1)) / d1;
        
        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat);
    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)
        % the thickness of unsaturated soil, (mm)
        d1 = zm(1);
        d2 = zgw - zm(1);
        
        % the unsaturated soil water, (mm)
        wa1_unsat = wa(1);
        wa2_unsat = (wa(2) * zm(2) - theta_sat * (zm(2) - d2)) / d2;
        
        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat);
    elseif zgw > zm(1) + zm(2)
        % the thickness of unsaturated soil, (mm)
        d1 = zm(1);
        d2 = zm(2);
        
        % the unsaturated soil water, (mm)
        wa1_unsat = wa(1);
        wa2_unsat = wa(2);
        
        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat);
    end

    if Pnet > 0.2 * Vmax
        srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax);
    else
        srf = 0;
    end
end

% actual water enter into the soil surface, (mm)
IWS = min(Vmax, Pnet - srf);
IWS = max(0, IWS);

% redundant water --> runoff (balance) # 限定边界条件
if IWS == Vmax || Vmax <= 0
    srf = Pnet - Vmax;
end

end
