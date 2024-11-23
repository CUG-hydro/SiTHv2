% Case 4 -- groundwater table in layer 4  %
function [wa, zgw, Tr, Es, uex] = swb_case4(wa, IWS, pEc, pEs, s_tem, s_vod, ...
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

    % unsaturated depth in layer #1~3
    d1 = zm(1);
    d2 = zm(2);
    d3 = zm(3);

    wa1 = wa(1);
    wa2 = wa(2);
    wa3 = wa(3);

    ks        = soilpar(1);  % hydraulic conductivity
    theta_sat = soilpar(3);  % saturated swc
    % theta_fc  = soilpar(5);  % field water capacity
    wwp       = soilpar(7);  % wilting point

    % ====== water supplement ====== %
    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1 % saturated
        wa1 = theta_sat;           % current soil water content
        vw1 = wc_s1 + IWS - wc_m1; % exceeded water
    else % un-saturated
        wa1 = wa1_unsat + IWS / d1; % soil water content in unsaturated zone
        vw1 = 0;
    end

    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = wa2;
    wc_s2 = d2 * wa2_unsat;

    % maximum water column in d2
    wc_m2 = d2 * theta_sat;

    if wc_s2 + vw1 >= wc_m2
        wa2 = theta_sat;           % current soil water content
        vw2 = wc_s2 + vw1 - wc_m2; % exceeded water to layer #3
    else
        wa2 = wa2_unsat + vw1 / d2; % soil water content in unsaturated zone
        vw2 = 0;                    % no exceeded water
    end

    % layer #3
    % existed water column in the unsaturated zone #3
    wa3_unsat = wa3;
    wc_s3 = d3 * wa3_unsat;

    % maximum water column in d3
    wc_m3 = d3 * theta_sat;

    if wc_s3 + vw2 >= wc_m3
        wa3 = theta_sat;           % current soil water content
        vw3 = wc_s3 + vw2 - wc_m3; % exceeded water
    else
        wa3 = wa3_unsat + vw2 / d3; % soil water content in unsaturated zone
        vw3 = 0;                    % no exceeded water
    end

    % ====== water consumption ====== %
    % Evapotranspiration %

    % distributed the potential T to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm);

    % Calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2, soilpar, pEc, pftpar);
    [f_sm3, ~] = swc_stress(wa3, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 .* s_vod * s_tem .* Tr_p1;
    Tr2 = f_sm2 .* s_vod * s_tem .* Tr_p2;
    Tr3 = f_sm3 .* s_vod * s_tem .* Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % Only consider about the first layer

    %% soil water drainage (unsaturated zone) %
    % ---------------------------------------------------------------- layer #1
    % update the soil moisture after ET, layer #1
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1 % wilting point
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % drainage from unsaturated zone, #1
    f1 = soil_drainage(wa1_unsat, theta_sat, ks, 0.048, 4.8);

    % update the soil moisture after drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);
    wa1 = min(wa1, 1); % < theta_s

    % ---------------------------------------------------------------- layer #2
    % update the soil moisture after ET, layer #2
    Tr2 = max(Tr2, 0); % reject negtive value
    Tr2 = min(Tr2, d2 * (wa2 - wwp)); % less than maximum avaliable water

    % drainage from unsaturated zone, #2
    f2 = soil_drainage(wa2_unsat, theta_sat, ks, 0.012, 1.2);

    % update the soil moisture after drainage, layer #2
    wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2;

    wa2 = max(wa2, 0); % > wwp 0
    wa2 = min(wa2, 1); % < theta_s 1

    if wa2 > theta_sat
        ff2 = (wa2 - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2 = theta_sat;
    else
        ff2 = 0;
    end

    % ---------------------------------------------------------------- layer #3
    % update the soil moisture after ET, layer #3, unsat-zone
    Tr3 = max(Tr3, 0); % reject negtive value
    Tr3 = min(Tr3, d3 * (wa3 - wwp)); % less than maximum avaliable water

    % drainage from unsaturated zone, #3
    f3 = soil_drainage(wa3_unsat, theta_sat, ks, 0.012, 1.2);

    % update the soil moisture after drainage, layer #3
    wa3 = (wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3;

    wa3 = max(wa3, 0); % > wwp
    wa3 = min(wa3, 1); % < theta_s
    if wa3 > theta_sat
        ff3 = (wa3 - theta_sat) * d3; % extra water from upper layer
        ff3 = max(ff3, 0);
        wa3_unsat = theta_sat;
        wa3 = theta_sat; % fix by kongdd, 2024-11-23
    else
        ff3 = 0;
        wa3_unsat = wa3;
    end

    %% The groundwater table depth %
    % total water recharge to groundwater
    F1 = f3 + ff3 + vw3;

    % total transpiration from groundwater
    Tr_g = 0;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes in groundwater table depth
    delta_zgw = delta_w / 0.2; % specific yield as 0.2
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % if F1 > 0
    %     fprintf("F1 = %.5f, R_sb = %.5f, delta_w = %.5f, delta_zgw = %.5f, zgw = %.5f \n", F1, R_sb, delta_w, delta_zgw, zgw);
    % end

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) && zgw <= zm(1) + zm(2) + zm(3)
        wa3 = (wa3_unsat * (zgw - zm(1) - zm(2)) + theta_sat * (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);
    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)
        wa2 = (wa2 * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;
    elseif zgw > 0 && zgw <= zm(1)
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
