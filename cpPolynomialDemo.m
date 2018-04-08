%% Demonstration of Polynomial
% Here we are plotting a variable specific heat as a function of
% temperature for air.

airCoeffs = fliplr([28.11, 0.1967e-2, 0.4802e-5,-1.966e-9]);
tmin = 273;
tmax = 1800;

temps = linspace(tmin, tmax);
cp = polyval(airCoeffs, temps);

plot(temps, cp);
xlabel('Temperature [k]');
ylabel('Specific Heat [kJ/kmol/k]');