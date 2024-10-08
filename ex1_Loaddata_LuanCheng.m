% Set spatial resolution
i = 522;
j = 2947;

% 0.1 deg
Soilraster = load('inpara\Soilraster.mat').Soilraster;
soil_type = Soilraster(i, j);

maskland = load('inpara\landmask01.mat').mask2;
maskland = maskland(i, j);

% Load the optimal temperature for plant growth
Topt = load('inpara\Topt.mat');
Topt = single(Topt.Topt_new);
Topt = Topt(i, j);

df = readtable("data/dat_栾城_ERA5L_1982-2019.csv");
nyear = 2016 - 1988 + 1;

%% model parameters 
Top = Topt;
PFTi = 22; % Land cover type
pftpar = get_pftpar(PFTi);
soilpar = get_soilpar(soil_type);

% retrieve data
dates = df.date;
years = year(dates);

%% 创建一个结构体，存放状态变量
% 对state敏感，KGE下降0.04, [0.43] -> [0.39]
state = struct();
uptval = load('data/State_spin-up_LuanCheng.mat').X_upti;
sm = uptval(1, 1:3); % 怀疑写颠倒了
zg = uptval(1, 4);
snowpack = uptval(1, 5);
state = mSiTH.update_state(state, sm, zg, snowpack);
