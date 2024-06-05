% Case 2 -- groundwater table in layer 2  %
function [wa, zgw, Tr, Es, uex] = swb_case2(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    %% INPUT:
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm

    % unsaturated depth in layer #1~2
    d1 = zm(1);
    d2 = zgw - d1;

    wa1 = wa(1);
    wa2 = wa(2);
    wa3 = wa(3);

    ks        = soilpar(1);  % hydraulic conductivity
    theta_sat = soilpar(3);  % saturated swc
    theta_fc  = soilpar(5);  % field water capacity
    wwp       = soilpar(7);  % wilting point

    % ====== water supplement ====== %
    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;
    wc_m1 = d1 * theta_sat; % maximum water column in d1

    if wc_s1 + IWS >= wc_m1
        wa1 = theta_sat;            % current soil water content  
        vw1 = wc_s1 + IWS - wc_m1;  % exceeded water
    else
        wa1 = wa1_unsat + IWS / d1; % soil water content in unsaturated zone
        vw1 = 0;                    % no exceeded water
    end
    
    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = (wa2 * zm(2) - theta_sat * (zm(2) - d2)) / d2;
    wc_s2 = d2 * wa2_unsat;
    wc_m2 = d2 * theta_sat; % maximum water column in d2

    if wc_s2 + vw1 >= wc_m2
        wa2_unsat = theta_sat;     % current soil water content
        vw2 = wc_s2 + vw1 - wc_m2; % exceeded water
    else
        wa2_unsat = wa2_unsat + vw1 / d2; % soil water content in unsaturated zone
        % calculate the adjusted swc#2 with considering the groundwater depth
        wa2 = (wa2_unsat * d2 + theta_sat * (zm(2) - d2)) / zm(2);
        vw2 = 0;
    end

    % layer #3
    % full filled with groundwater

    % ====== water consumption ====== %
    %% Evapotranspiration %
    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm);

    % divide Tr_p2 into unsaturated zone and saturated zone
    Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + ...
    (zm(2) - d2) * theta_sat);
    Tr_p2_g = Tr_p2 * ((zm(2) - d2) * theta_sat) / (d2 * wa2_unsat + ...
        (zm(2) - d2) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2_unsat, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 * s_vod * s_tem * Tr_p1;
    Tr2_u = f_sm2 * s_vod * s_tem * Tr_p2_u;
    Tr2_g = s_vod * s_tem * Tr_p2_g;
    Tr2 = Tr2_u + Tr2_g;
    Tr3 = s_vod * s_tem * Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only considering the first layer

    %% soil water drainage (unsaturated zone) %
    % layer #1
    % check Es & Tr, layer #1   unsat-zone
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1  % wilting point (0) for first layer
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % Ducharne % Polcher, 1998
    dd = 1.5;
    Dmin = 0.048; % mm day-1
    Dmax = 4.8; % mm day-1

    thx1 = wa1 / theta_sat;

    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end

    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);

    % update the soil moisture after drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);

    % ---------------------------------------------------------------- layer #2
    % check Tr, layer #2   unsat-zone
    Tr2_u = max(Tr2_u, 0, d2 * (wa2_unsat - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx2 = wa2_unsat / theta_sat;

    if thx2 < 0.75
        Perc2 = Dmin * thx2;
    else
        Perc2 = Dmin * thx2 + (Dmax - Dmin) * thx2^dd;
    end

    % drainage from unsaturated zone, #2
    f2 = min(ks, Perc2);

    % update the soil moisture after ET & drainage, layer #2 un
    wa2_unsat = (wa2_unsat * d2 + f1 - f2 - Tr2_u) / d2;
    wa2_unsat = clamp(wa2_unsat, 0, 1); % > 0

    if wa2_unsat > theta_sat
        ff2 = (wa2_unsat - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2_unsat = theta_sat;
    else
        ff2 = 0;
    end

    % layer #3
    % full filled with groundwater

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f2 + ff2 + vw2;

    % total transpiration from groundwater
    Tr_g = Tr2_g + Tr3;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
%     if delta_w < 0 % decline of the water table
% 
%         delta_zgw = delta_w / (theta_sat - theta_fc);
% 
%     else % increase of the water table
    delta_zgw = delta_w / (theta_sat - (wa1 + wa2_unsat) / 2);
%     end

    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)
        wa2 = (wa2_unsat * d2 + theta_fc * (zm(2) - d2)) / zm(2);
        wa3 = theta_fc;
    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)
        wa2 = (wa2_unsat * d2 + theta_fc * (zm(2) - d2)) / zm(2);
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);
    elseif zgw > zm(1) && zgw < zm(1) + zm(2)
        wa2 = (wa2_unsat * (zgw - zm(1)) + theta_sat * ...
            (zm(1) + zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;
    elseif zgw > 0 && zgw < zm(1)
        wa1 = (wa1 * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;
    elseif zgw <= 0
        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;
        uex = -zgw * theta_fc; % excess water to soil surface, mm
    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);
end
