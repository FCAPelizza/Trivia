clear all
% syms kg m K s
% mu = kg/(m*s);
% kt = kg*m/(s^3*K);
% D = m^2/s;
% P = kg/m/s^2;
% time = mu/P
% temperature = P*D/kt
% length = (D*mu/P)^0.5
% mass = mu*(D*mu/P)^0.5*mu/P
% speed = length/time % (D*P/mu)^0.5/(mu/P)
%syms mu kt D P
P = 11034*9.81*1000; % pressure at the bottom of the mariana trench
mu = 10^19; % estimation of viscosity of glass at room temperature
kt = 0.00565; % thermal conductivity of Xenon
D = 1.3*10^-30; % diffusivity of Al in Cu
time = mu/P % around 3100 years
temperature = P*D/kt %
length = (D*mu/P)^0.5 % around 0.1 nanometer, the diameter of atomic helium
mass = mu*(D*mu/P)^0.5*mu/P % around the mass of pluto
speed = length/time % it would take you 10'000 ages of the universe to travel a meter
acceleration = speed/time
force = acceleration * mass % this is 10^7 less force than produced by a snail in acceleration
energy = force * length % 10^5 less energy than the kinetic energy of a dust particle travelling at 0.01 m/s
power = energy / time
%pretty(power)
%pretty(speed)