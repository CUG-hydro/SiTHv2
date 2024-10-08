%% 栾城站: 模型测试
% kongdd, CUG, 2024-10-08
clc, clear;
ex1_LuanCheng;

N = size(df, 1);
ET = zeros(N, 1);
Tr = zeros(N, 1);
Es = zeros(N, 1);
Ei = zeros(N, 1);
Esb = zeros(N, 1); % Snow sublimation
RF = zeros(N, 1);
GW = zeros(N, 1);
R = zeros(N, 1);
SM = zeros(N, 3);

for yr = 1989:2016
    k = yr - 1988 + 1; % 因为1988是第一年，索引从1开始
    spinfg = 0;
    fprintf('[%d, spinfg=0] ... \n', yr)
    
    inds = find(years == yr);
    days = length(inds); % 获取该年的天数
    
    d = df(inds, :);
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
    
    [ET(inds), Tr(inds), Es(inds), Ei(inds), Esb(inds), SM(inds, :), RF(inds), GW(inds), state] = ...
        cal_SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, Gi, LAI, Top, s_VODi, ...
        soilpar, pftpar, state, spinfg);
end

SM1 = SM(:,1);
SM2 = SM(:,2);
SM3 = SM(:,3);
df_out = table(dates, ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3);
writetable(df_out, "data/OUTPUT_栾城.csv")
