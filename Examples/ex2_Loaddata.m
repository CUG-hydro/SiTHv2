% 针对CRO进行参数优化
Soilraster = load('inpara\Soilraster.mat').Soilraster;


% three global variables
st = readtable("data/siteInfo_CRO_5sp.csv");
ra_soil = load('inpara\Soilraster.mat').Soilraster; % LAT, LON
ra_Topt = single(load('inpara\Topt.mat').Topt_new);

k = 1;
[soilpar, Topt, pftpar, state] = init_param(st, ra_soil, ra_Topt, k);


% 根据土壤属性，对参数进行初始化。
% 一开始对田间进行漫灌，各层从fc开始。
function [soilpar, Topt, pftpar, state] = init_param(st, ra_soil, ra_Topt, k)
    % global LON LAT ra_soil ra_Topt
    cellsize=0.1;
    LON = -180+cellsize/2:cellsize:180;
    LAT = flip(-90+cellsize/2:cellsize:90);

    x = st.lon(k);
    y = st.lat(k);
    [i, j] = findnear(LON, LAT, x, y);

    soil_type = ra_soil(j, i);
    soilpar = get_soilpar(soil_type);

    Topt = ra_Topt(j, i);

    PFTi = 22; % Land cover type
    pftpar = get_pftpar(PFTi);

    theta_sat = soilpar(3); % saturated wa
    sm = ones(3, 1) * theta_sat;
    zg = 0;
    snowpack = 0;
    state = mSiTH.update_state(struct(), sm, zg, snowpack);
end

function [i, j] = findnear(LON, LAT, x, y)
    [~, i] = min(abs(LON - x));
    [~, j] = min(abs(LAT - y));
end
