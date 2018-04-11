function [] = plotTS(T,P,title)
    delta_s = findDelta(T,P);
    delta_s = [delta_s; 0];
    plot(delta_s,[T T(1)],'b*--')
    box off
    ylabel('Temperature (K)')
    xlabel('Entropy (kJ/kg*K)')
    title(title)
end