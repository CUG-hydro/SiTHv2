% Case 3 -- groundwater table in layer 3  %
function [wa, zgw, Tr, Es, uex] = swb_case3(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % unsaturated depth in layer #1
    d1 = zm(1);
    % unsaturated depth in layer #2
    d2 = zm(2);
    % unsaturated depth in layer #3
    d3 = zgw - zm(1) - zm(2);

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1

        % current soil water content
        wa1 = theta_sat;
        % exceeded water
        vw1 = wc_s1 + IWS - wc_m1;
    else

        % soil water content in unsaturated zone
        wa1 = wa1_unsat + IWS / d1;
        vw1 = 0; % no exceeded water
    end
   
    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = wa2;
    wc_s2 = d2 * wa2_unsat;

    % maximum water column in d2
    wc_m2 = d2 * theta_sat;

    if wc_s2 + vw1 >= wc_m2

        % current soil water content
        wa2 = theta_sat;
        % exceeded water
        vw2 = wc_s2 + vw1 - wc_m2;
    else
        % soil water content in unsaturated zone
        wa2 = wa2_unsat + vw1 / d2;
        vw2 = 0; % no exceeded water
    end

    % layer #3
    % existed water column in the unsaturated zone #3
    wa3_unsat = (wa3 * zm(3) - theta_sat * (zm(3) - d3)) / d3;
    wc_s3 = d3 * wa3_unsat;

    % maximum water column in d3
    wc_m3 = d3 * theta_sat;

    if wc_s3 + vw2 >= wc_m3
        wa3 = theta_sat;
        vw3 = wc_s3 + vw2 - wc_m3;
    else
        % soil water content in unsaturated zone
        wa3_unsat = wa3_unsat + vw2 / d3;
        % calculate the adjusted swc#3 with considering the groundwater depth
        wa3 = (wa3_unsat * d3 + theta_sat * (zm(3) - d3)) / zm(3);
        % no exceeded water
        vw3 = 0;
    end

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % divide Tr_p3 into unsaturated and saturated zone at the layer #3
    Tr_p3_u = Tr_p3 * (d3 * wa3_unsat) / (d3 * wa3_unsat + ...
    (zm(3) - d3) * theta_sat);
    Tr_p3_g = Tr_p3 * ((zm(3) - d3) * theta_sat) / (d3 * wa3_unsat + ...
        (zm(3) - d3) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2, soilpar, pEc, pftpar);
    [f_sm3, ~] = swc_stress(wa3_unsat, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 * s_vod * s_tem * Tr_p1;
    Tr2 = f_sm2 * s_vod * s_tem * Tr_p2;
    Tr3_u = f_sm3 * s_vod * s_tem * Tr_p3_u;
    Tr3_g = s_vod * s_tem * Tr_p3_g;
    Tr3 = Tr3_u + Tr3_g;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only considering the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %
    % ---------------------------------------------------------------- layer #1
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1 % wilting point -- 0, first layer
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % Ducharne & Polcher, 1998
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

    % update the soil moisture after ET & drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);

    % ---------------------------------------------------------------- layer #2
    Tr2 = max(Tr2, 0); % reject negtive value
    Tr2 = min(Tr2, d2 * (wa2 - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx2 = wa2 / theta_sat;

    if thx2 < 0.75
        Perc2 = Dmin * thx2;
    else
        Perc2 = Dmin * thx2 + (Dmax - Dmin) * thx2^dd;
    end

    % drainage from unsaturated zone, #2
    f2 = min(ks, Perc2);

    % update the soil moisture after ET & drainage, layer #2
    wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2;

    wa2 = max(wa2, 0); % > wwp--0
    wa2 = min(wa2, 1); % < 1

    if wa2 > theta_sat
        ff2 = (wa2 - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2 = theta_sat;
    else
        ff2 = 0;
    end

    % ---------------------------------------------------------------- layer #3
    % check Tr, layer #3   unsat-zone
    Tr3_u = max(Tr3_u, 0); % reject negtive value
    Tr3_u = min(Tr3_u, d3 * (wa3_unsat - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx3 = wa3_unsat / theta_sat;

    if thx3 < 0.75
        Perc3 = Dmin * thx3;
    else
        Perc3 = Dmin * thx3 + (Dmax - Dmin) * thx3^dd;
    end

    % drainage from unsaturated zone, #3
    f3 = min(ks, Perc3);

    % update the soil moisture after ET & drainage, layer #3
    wa3_unsat = (wa3_unsat * d3 + f2 + ff2 - f3 - Tr3_u) / d3;

    wa3_unsat = max(wa3_unsat, 0); % > wwp 0
    wa3_unsat = min(wa3_unsat, 1); % < theta_s 1
    
    if wa3_unsat > theta_sat
        ff3 = (wa3_unsat - theta_sat) * d3; % extra water from upper layer
        ff3 = max(ff3, 0);
        wa3_unsat = theta_sat;
    else
        ff3 = 0;
    end

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f3 + ff3 + vw3;

    % total transpiration from groundwater
    Tr_g = Tr3_g;

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
    delta_zgw = delta_w / (theta_sat - (wa1 + wa2 + wa3_unsat) / 3);
%     end
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)
        wa3 = (wa3_unsat * d3 + theta_fc * (zm(3) - d3)) / zm(3);
    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)
        wa3 = (wa3_unsat * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);
    elseif zgw > zm(1) && zgw < zm(1) + zm(2)
        wa2 = (wa2 * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) ...
            / zm(2);
        wa3 = theta_sat;
    elseif zgw > 0 && zgw < zm(1)
        wa1 = (wa1 * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;
    elseif zgw <= 0
        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;
        uex = -zgw * theta_sat; % excess water to soil surface, mm
    end
    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);
end
