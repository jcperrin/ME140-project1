function [s] = findDelta(T,P)
    % Finds difference in entropy given temperature, pressure and variable cp
    if length(T) ~= length(P)
        ME = MException('T and P are not the same length. Please retry');
        throw(ME)
    end

    % Constants
    R = .287;
    mm_air = 28.97; %kg/kmol
    airCoeffs = fliplr([28.11, 0.1967e-2, 0.4802e-5,-1.966e-9]);

    s = zeros(length(T),1);

    % Find c_p for given state
    for i= 1:length(T)-1
        T1 = T(i);
        T2 = T(i+1);

        P1 = P(i);
        P2 = P(i+1);
        scircley = airCoeffs(4)*log(T2/T1) + airCoeffs(3)*(T2-T1) + ...
            .5*airCoeffs(2)*(T2^2 - T1^2) + (1/3)*airCoeffs(1)*(T2^3-T1^3);
        s(i+1) = s(i) + scircley/mm_air - R*log(P2/P1);
    end

end