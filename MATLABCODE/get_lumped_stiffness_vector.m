function [EIs] = get_lumped_stiffness_vector(x)
    % x is a vector with the position of each station
    % EIs is the vector with corresponding lumped bending stifness for each station

    % Position of stations [m]
    x_of_stations = [0; 18; 36; 54; 108; 144; 157.5];

    % Bending Stiffness distribution [N.m^2]
    EI_data = (10^9) * [3543; 5905; 7086; 7676; 7676; 6495; ...
                        4724];

    % Regression to obtain the polynomial responsible for describing
    % the bending stiffness at any x
    EI_fit = fit(x_of_stations, EI_data, 'linearinterp');

    EIs = zeros(length(x), 1);

    for i = 1:(length(x))
        EIs(i) = EI_fit(x(i));
    end

end
