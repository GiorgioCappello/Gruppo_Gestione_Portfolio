function [PI_new, sigma_new, w_new] = calculate_view(sigma_SCC, PI, P, Q, Certainty)
    % Calcolo di una view generica
    % Calcolo di Omega
    for i=1:size(P,1)
        Omega(i,i)=P(i,:)*sigma_SCC*P(i,:)'*((1-Certainty(i))*1/Certainty(i));
    end
    
    % Calcolo dei nuovi parametri
    PI_new = inv(inv(sigma_SCC) + P' * inv(Omega) * P) * (inv(sigma_SCC) * PI + P' * inv(Omega) * Q);
    sigma_new = inv(inv(sigma_SCC) + P' * inv(Omega) * P);
    w_new = (sigma_new \ PI_new) / sum(sigma_new \ PI_new);
end
