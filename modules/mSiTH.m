classdef mSiTH
    %SiTH model class
    properties
        
    end
    
    methods (Static)
        function [State] = update_state(State, sm, zg, snowpack)
            % Create a structure to store state variables
            % 
            % 状态变量需要连续，传递到下一年中；模型对初始状态state敏感。
            % - sm: 采用warming-up的方式获取，warming-up period可取3年
            % - zg: 采用spin-up的方式获取，spin-up period可取100年
            % 
            %% Argument Specification
            % - sm: soil water content in three layers, [m^3 m^-3]
            % - zg: groundwater depth, [mm]
            % - snowpack: snowpack depth
            sm(sm < 0) = 0.01; % set the minimum value for soil moisture
            State.sm = sm;
            State.zg = zg;
            State.snowpack = snowpack;
        end
        
        function [sm, zg, snowpack] = get_state(State)
            % Get the state variables
            %% Argument Specification
            % - State: a structure to store state variables
            sm = State.sm;
            zg = State.zg;
            snowpack = State.snowpack;
        end
    end
end
