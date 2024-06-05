function [x] = clamp(x, min, max)
    if x < min
        x = min;
    elseif x > max
        x = max;
    end
end
