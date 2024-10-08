% run SiTHv2 model
function [ETs, Trs, Ess, Eis, Esbs, SM, RF, GW, state] = cal_SiTHv2_site(Rni, Tai, Tasi, Precii,...
    Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, spinfg)
%% INPUT
% - `Top`   : optimal growth temperature for plant, degC
% - `wa`    : Soil moisture (last step)
% - `zgw`   : groundwater table depth, mm
% - `snp`   : Snow package (old), mm day-1
% - `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
[wa, zgw, snp] = mSiTH.get_state(state);
ntime = size(Rni, 1);

if spinfg == 1 % spin-up
    ETs = zeros(ntime, 1);
    Trs = zeros(ntime, 1);
    Ess = zeros(ntime, 1);
    Eis = zeros(ntime, 1);
    Esbs = zeros(ntime, 1);
    SM = zeros(ntime, 3);
    RF = zeros(ntime, 1);
    GW = zeros(ntime, 1);
    
    for k = 1 : 100 % set the spin-up time (100 years)
        for i = 1:ntime
            Rn = Rni(i, 1);
            Ta = Tai(i, 1);
            Tas = Tasi(i, 1);
            Pe = Precii(i, 1);
            Pa = Pai(i, 1);
            G = Gi(i, 1);
            LAI = LAIii(i, 1);
            s_VOD = s_VODi(i, 1);
            
            [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, ~, ~, ~] = SiTHv2(Rn, Ta, Tas, Top, ...
                Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp);
            
            ETs(i, 1) = Et;
            Trs(i, 1) = Tr;
            Ess(i, 1) = Es;
            Eis(i, 1) = Ei;
            Esbs(i, 1) = Esb;
            SM(i, :) = wa;
            RF(i, 1) = srf;
            GW(i, 1) = zgw;
        end
    end
else 
    % normal calculation
    ETs = zeros(ntime, 1);
    Trs = zeros(ntime, 1);
    Ess = zeros(ntime, 1);
    Eis = zeros(ntime, 1);
    Esbs = zeros(ntime, 1);
    SM = zeros(ntime, 3);
    RF = zeros(ntime, 1);
    GW = zeros(ntime, 1);
    
    for i = 1 : ntime
        Rn = Rni(i, 1);
        Ta = Tai(i, 1);
        Tas = Tasi(i, 1);
        Pe = Precii(i, 1);
        Pa = Pai(i, 1);
        G = Gi(i, 1);
        LAI = LAIii(i, 1);
        s_VOD = s_VODi(i, 1);
        
        [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, ~, ~, ~] = SiTHv2(Rn, Ta, Tas, Top, ...
            Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp);
        
        ETs(i, 1) = Et;
        Trs(i, 1) = Tr;
        Ess(i, 1) = Es;
        Eis(i, 1) = Ei;
        Esbs(i, 1) = Esb;
        SM(i, :) = wa;
        RF(i, 1) = srf;
        GW(i, 1) = zgw;
    end
end
state = mSiTH.update_state(state, wa, zgw, snp);

end
