function [m] = get_lumped_mass_vector(x)
    % x is a vector with the position of each station
    % mi is the vector with corresponding lumped masses for each station
    % we consider the added mass here

    % Position of stations
    x_of_stations = [0; 7.88; 15.75; 23.63; 31.50; 39.38; 47.25; 55.13; 63.00; 70.88; ...
                    78.75; 86.83; 94.50; 102.38; 110.25; 118.13; 126.00; 133.88; 141.75; ...
                    149.63; 157.50];

    % Mass distribution [kg/m]
    % From table on project assignment
    % W_table_data = 1000 * [17; 38; 59; 73; 88; 103; 117; 132; 146; 146; 146; 146; 146; ...
    %                         135; 123; 111; 99; 88; 76; 51; 27];
    W_table_data = 1000 * [17; 38; 59; 73; 88; 103; 117; 132; 146; 146; 146; 146; 146; ...
                            135; 123; 111; 99; 88; 76; 51; 27];

    % Result with ANSYS for j=2
    % Ideally, we should use this:
    % get_added_mass = @(x) (pi / 8) * p * B(x) * C(x) * Jn;
    % and consider the right B(x), C(x) and Jn for each case
    W_added_mass_data = 1000 * [0.91; 28.35; 57.76; 91.25; 124.30; 159.33; 186.61; 198.36; ...
                                211.87; 222.42; 222.42; 203.31; 172.46; 145; 107.19; ...
                                67.15; 37.70; 17.61; 5.75; 2.00; 0];

    W_total_data = W_table_data + W_added_mass_data;

    % Regression to obtain the polynomial responsible for describing weigth per length against x
    W = fit(x_of_stations, W_table_data, 'linearinterp'); % Weight per length [kg/m]

    m = zeros(length(x), 1);

    for i = 1:(length(x) - 1)
        x_left = x(i); x_right = x(i + 1);
        dx = x_right - x_left;

        xx = [x_left, x_right];
        M = trapz(xx, [W(x_left), W(x_right)]);
        x_CG = (dx / 3) * ((2 * W(x_right) + W(x_left)) / (W(x_right) + W(x_left)));

        m_left = M * (x_right + 2 * x_left - x_CG) / (x_right + x_left);
        m_right = M * (x_CG - x_left) / (x_right + x_left);

        m(i) = m(i) + m_left;
        m(i + 1) = m_right;
    end

end
