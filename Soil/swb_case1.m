% Case 1 -- groundwater table in layer 1  %
function [wa, zgw, Tr, Es, uex] = swb_case1(wa, IWS, pEc, pEs, s_tem, s_vod, ...
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

    % unsaturated depth in layer #1
    d1 = zgw;

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
    wa1_unsat = (wa1 * zm(1) - theta_sat * (zm(1) - d1)) / d1; % 未饱和区域SM, [m3 m-3]
    wc_s1 = d1 * wa1_unsat; % [mm]
    wc_m1 = d1 * theta_sat; % maximum water column in d1, [mm]

    if wc_s1 + IWS >= wc_m1
        wa1_unsat = theta_sat;      % current soil water content
        vw1 = wc_s1 + IWS - wc_m1;  % exceeded water
    else
        % soil water content in unsaturated zone
        wa1_unsat = wa1_unsat + IWS / d1; 
        % calculate the adjusted swc#1 with considering the groundwater depth
        wa1 = (wa1_unsat * d1 + theta_sat * (zm(1) - d1)) / zm(1);
        vw1 = 0; 
    end

    % layer #2 - saturated, full filled with groundwater
    % layer #3 - saturated

    % ====== water consumption ====== %
    % Evapotranspiration %
    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm);

    % divide Tr_p1 into unsaturated zone and saturated zone
    Tr_p1_u = Tr_p1 * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm(1) - d1) * theta_sat);
    Tr_p1_g = Tr_p1 * ((zm(1) - d1) * theta_sat) / (d1 * wa1_unsat + (zm(1) - d1) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1_u = f_sm1 * s_vod * s_tem * Tr_p1_u;
    Tr1_u = clamp(Tr1_u, 0, d1 * (wa1_unsat - wwp)); % less than maximum avaliable water
    
    Tr1_g = s_vod * s_tem * Tr_p1_g;
    Tr1 = Tr1_u + Tr1_g;
    Tr2 = s_vod * s_tem * Tr_p2;
    Tr3 = s_vod * s_tem * Tr_p3;

    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only consider about the first layer
    Es_u = Es * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm(1) - d1) * theta_sat);
    Es_u = clamp(Es_u, 0, d1 * wa1_unsat);

    %% soil water drainage (unsaturated zone) %
    % layer #1
    % drainage from unsaturated zone, #1
    f1 = soil_drainage(wa1_unsat, theta_sat, ks, 0.048, 4.8);

    % update the soil moisture after ET & drainage, layer #1
    wa1_unsat = (wa1_unsat * d1 - f1 - Es_u - Tr1_u) / d1;
    wa1_unsat = max(wa1_unsat, 0);

    % layer #2   full filled with groundwater
    % layer #3   full filled with groundwater

    %% The groundwater table depth
    % total water recharge to groundwater
    F1 = f1 + vw1;

    % total transpiration from groundwater
    Tr_g = Tr1_g + Tr2 + Tr3;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
%     if delta_w < 0 % decline of the water table
%         delta_zgw = delta_w / theta_fc;
%     else % increase of the water table
      delta_zgw = delta_w / (theta_sat - wa1_unsat);
%     end
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)
        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = theta_fc;
    elseif zgw > zm(1) + zm(2) && zgw <= zm(1) + zm(2) + zm(3)
        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);
    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)
        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = (theta_fc * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;
    elseif zgw > 0 && zgw <= zm(1)
        wa1 = (wa1_unsat * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
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
