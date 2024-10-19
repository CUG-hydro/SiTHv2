%% 栾城站: 模型测试
% kongdd, CUG, 2024-10-08
clc, clear;

d = readtable("data/dat_栾城_ERA5L_2010.csv");
spinfg = 0;

Rn   = d.Rn;
Ta   = d.Tavg;
Prcp = d.Prcp;
Pa   = d.Pa;
LAI  = d.LAI;
VOD  = d.VOD;

Tas = Ta;            % 有效积温
Tas(Tas < 0) = 0;    % 去除小于0的值
Tas = cumsum(Tas);

Gi = 0.4 .* Rn .* exp(-0.5 .* LAI); % G_soil
s_VODi = (VOD ./ max(VOD)).^0.5; % VOD-stress


zgws = [0, 25, 1000, 2000, 6000.0];
for i = 4:length(zgws)
    zgw = zgws(i);
    [soilpar, Top, pftpar, state] = init_param();
    state.zg = zgw;

    [ET, Tr, Es, Ei, Esb, SM, RF, GW, state] = SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, Gi, LAI, Top, s_VODi, ...
    soilpar, pftpar, state, spinfg);

    SM1 = SM(:,1);
    SM2 = SM(:,2);
    SM3 = SM(:,3);
    df_out = table(ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3);
    
    fout = sprintf("./Examples/OUTPUT/OUTPUT_栾城_2010_MATLAB_zgw=%d.csv", zgw);
    writetable(df_out, fout)
end
% 28*10年
% 前10次spin-up，最后一次取结果

% system("Rscript data/Figure1_栾城_SiTH.R");

% 根据土壤属性，对参数进行初始化。
% 一开始对田间进行漫灌，各层从fc开始。
function [soilpar, Topt, pftpar, state] = init_param()
    soil_type = 2;
    soilpar = get_soilpar(soil_type);

    Topt = 24.0;

    PFTi = 22; % Land cover type
    pftpar = get_pftpar(PFTi);

    theta_sat = soilpar(3); % saturated wa
    sm = ones(3, 1) * theta_sat;
    zgw = 0;
    snowpack = 0;
    state = mSiTH.update_state(struct(), sm, zgw, snowpack);
end
