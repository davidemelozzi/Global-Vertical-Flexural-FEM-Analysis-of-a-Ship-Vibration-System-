function [Y] = extract_mode_from_station_vector(v)
    Y = zeros(length(v), 1);

    for i = 1:length(v)
        Y(i) = v{i}(1);
    end

end
