function [Y, wn] = analytical(x, EI, p, A, L)
    beta1_times_L = 4.7300041;
    beta2_times_L = 7.853205;
    beta3_times_L = 10.995608;
    beta1 = (beta1_times_L / L);
    beta2 = (beta2_times_L / L);
    beta3 = (beta3_times_L / L);

    %  Calculating first three natural frequencies
    square_root = sqrt(EI / (p * A * (L^4)));
    wn = [(beta1_times_L^2) * square_root,
        (beta2_times_L^2) * square_root,
        (beta3_times_L^2) * square_root];

    %  Calculating first three eigenfucntions/mode shapes
    function [Y] = mode_shape(x, beta, L)
        Cn = 1;
        alfa_n = (sinh(beta * L) - sin(beta * L)) / (cos(beta * L) - cosh(beta * L));
        Y = Cn * (sinh(beta * x) - sin(beta * x)) + alfa_n * (cosh(beta * x) - cos(beta * x));
    end

    Y = [mode_shape(x, beta1, L),
        mode_shape(x, beta2, L),
        mode_shape(x, beta3, L)];
end
