function [State] = update_state(State, SM, ZG, snowpack)
  % Create a structure to store state variables
  %% Argument Specification
  % - SM: soil water content in three layers, [m^3 m^-3]
  % - ZG: groundwater depth, [mm]
  % - snowpack: snowpack depth
  SM(SM < 0) = 0.01; % set the minimum value for soil moisture
  State.SM = SM;
  State.ZG = ZG;
  State.snowpack = snowpack;
end
