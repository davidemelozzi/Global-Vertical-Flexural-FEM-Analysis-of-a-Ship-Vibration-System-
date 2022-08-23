function [wn] = kumai(B, D, T, Lpp, displacement)
    % Calculate the natural frequencies by using the simple Kumai's formula
    % B = Beam
    % D = Depth
    % T = Draft
    % Lpp = Length between perpendiculars
    % displacement = displacement of ship

    % Area moment of inertia of ship's midsection [m^4]
    % Depth Kumai correction
    D = D * 0.5;
    I = B * (D^3) / 12;

    % Displacement including virtual mass [t]
    displacement_with_vmass = (1.2 + B / (3 * T)) * displacement;

    % First natural frequency [hz]
    wn(1) = (1/60) * (3.07 * 10^6) * sqrt(I / (displacement_with_vmass * Lpp^3));

    num_of_natural_modes = 3;
    alpha = 1.02; % considering a tankero ship

    for n = 3:(num_of_natural_modes + 1)
        wn(n - 1) = wn(1) * (n - 1)^alpha;
    end

end
