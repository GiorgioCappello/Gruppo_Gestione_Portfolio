function [PI_new, sigma_new, w_new] = calculateView(sigma_SCC, PI, P, Q, Certainty)
    % Calcolo di Omega
    Omega = P * sigma_SCC * P' * ((1 - Certainty) / Certainty);
    
    % Calcolo dei nuovi parametri
    PI_new = inv(inv(sigma_SCC) + P' * inv(Omega) * P) * (inv(sigma_SCC) * PI + P' * inv(Omega) * Q);
    sigma_new = inv(inv(sigma_SCC) + P' * inv(Omega) * P);
    w_new = (sigma_new \ PI_new) / sum(sigma_new \ PI_new);
end
