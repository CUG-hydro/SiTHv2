% Case 0 -- groundwater overflow %
function [wa, zgw, Tr, Es, uex] = swb_case0(wa, IWS, pEc, pEs, s_tem, s_vod, ...
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

    % all saturated in layer #1 #2 #3
    % old soil water content in layer 1-3
    wa1 = wa(1);
    wa2 = wa(2);
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    % ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    % wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1 - saturated
    % full filled with groundwater

    % layer #2 - saturated
    % full filled with groundwater

    % layer #3 - saturated
    % full filled with groundwater

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % actual transpiration
    Tr1 = s_vod * s_tem * Tr_p1;
    Tr2 = s_vod * s_tem * Tr_p2;
    Tr3 = s_vod * s_tem * Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;
    % actual soil evaporation
    Es = pEs; % only consider about the first layer

    % layer #1   full filled with groundwater
    % layer #2   full filled with groundwater
    % layer #3   full filled with groundwater

    %% The groundwater table depth %
    F1 = IWS;  % total water recharge to groundwater from unsaturated zone
    Tr_g = Tr; % total transpiration from groundwater
    
    R_sb_max = 39; % groundwater discaharge, [mm day-1]
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb; % should below zero

    % changes of groundwater table depth
    delta_zgw = delta_w / (theta_sat - theta_fc); % 按照田间持水量进行排水
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)
        wa1 = theta_fc;
        wa2 = theta_fc;
        wa3 = theta_fc;
    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)
        wa1 = theta_fc;
        wa2 = theta_fc;
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3); % 水面之下饱和、落水之上fc
    elseif zgw > zm(1) && zgw < zm(1) + zm(2)
        wa1 = theta_fc;
        wa2 = (theta_fc * (zgw - zm(1)) + theta_sat * (zm(1) + ...
            zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;
    elseif zgw > 0 && zgw < zm(1)
        wa1 = (theta_fc * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;
    elseif zgw <= 0
        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;
        uex = -zgw * theta_fc; % excess water to soil surface, mm
        % uex = -zgw * (theta_sat - theta_fc); % excess water to soil surface, mm
    end
    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);
end
