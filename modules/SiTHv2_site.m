% run SiTHv2 model
function [ETs, Trs, Ess, Eis, Esbs, SM, RF, GW, state] = SiTHv2_site(Rni, Tai, Tasi, Precii,...
    Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, spinfg)
%% INPUT
% - `Top`   : optimal growth temperature for plant, degC
% - `state`: 
%   + `wa`    : Soil moisture (last step)
%   + `zgw`   : groundwater table depth, mm
%   + `snp`   : Snow package (old), mm day-1
% - `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
ntime = size(Rni, 1);

output = struct();
output.ETs  = zeros(ntime, 1);
output.Trs  = zeros(ntime, 1);
output.Ess  = zeros(ntime, 1);
output.Eis  = zeros(ntime, 1);
output.Esbs = zeros(ntime, 1);
output.SM   = zeros(ntime, 3);
output.RF   = zeros(ntime, 1);
output.GW   = zeros(ntime, 1);

if spinfg == 1 % spin-up
    for k = 1 : 100 % set the spin-up time (100 years)
        [state, output] = run_model(Rni, Tai, Tasi, Precii,...
            Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output);
    end
else
    [state, output] = run_model(Rni, Tai, Tasi, Precii,...
            Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output);
end

ETs  = output.ETs;
Trs  = output.Trs;
Ess  = output.Ess;
Eis  = output.Eis;
Esbs = output.Esbs;
SM   = output.SM;
RF   = output.RF;
GW   = output.GW;
end


function [state, output] = run_model(Rni, Tai, Tasi, Precii, ...
    Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output)

[wa, zgw, snp] = mSiTH.get_state(state);
ntime = size(Rni, 1);

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
    
    output.ETs(i, 1) = Et;
    output.Trs(i, 1) = Tr;
    output.Ess(i, 1) = Es;
    output.Eis(i, 1) = Ei;
    output.Esbs(i, 1) = Esb;
    output.SM(i, :) = wa;
    output.RF(i, 1) = srf;
    output.GW(i, 1) = zgw;
end
state = mSiTH.update_state(state, wa, zgw, snp);
end
