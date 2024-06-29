clear
clc
close all
%% main data

CeramicMug_ThermalConductivity = 3.8; %W/(m*K)
Steam_DynamicViscosity = 0.0125*10^-3; %Pa*s
Steam_Diffusivity = 0.325*10^-4; %m^2/s
Air_ThermalConductivity = 26.24e-03; % W/(m*K)
Water_MassDensity = 1000; %kg/m^3
Steam_SpecificHeat = 1000; %J/(kg*K)
Mug_Thickness = 5e-03; %m
Water_SpecificHeat = 4186; %J/(kg*K)
Water_HeatOfVaporisation = 2300e+03; %J/kg
Air_MassDensity = 1.184; %kg/m^3
Pressure = 101325; %Pa
Water_MolarMass = 18e-03; %kg/mol
Air_MolarMass = 28.96e-03; %kg/mol
MugDiameter = 0.04;
Air_RelativeHumidityPercentage = 30;
Air_RelativeHumidity = Air_RelativeHumidityPercentage/100;
Air_TemperatureCelsius = 20;
Air_Temperature = Air_TemperatureCelsius + 273.15; %K

%% initial Conditions

IntialTemperatureCelsius = 60; %C
InitialLiquidHeight = 0.001; %m
TimeEndIntegrationMin = 3000; % min

InitialMass = InitialLiquidHeight*pi*MugDiameter^2/4*Water_MassDensity;
odeoptions = odeset('NonNegative',1);

IntialConditions = [InitialMass,IntialTemperatureCelsius+273.15];

[t,Ys] = ode89(@(t,inputs)odefun(t,inputs,Water_MassDensity,MugDiameter,Steam_DynamicViscosity...
    ,Steam_SpecificHeat,Air_ThermalConductivity,Steam_Diffusivity...
    ,CeramicMug_ThermalConductivity,Water_HeatOfVaporisation,Water_SpecificHeat,Air_RelativeHumidity...
    ,Air_Temperature,Pressure,Water_MolarMass,Air_MolarMass,Mug_Thickness),[0,TimeEndIntegrationMin*60],IntialConditions,odeoptions);



figure(1)

plot(t/60,Ys(:,1)/(pi*MugDiameter^2/4*Water_MassDensity)*100);
title('Height of the liquid in time')
ylabel('Height of the liquid [cm]')
title('time [min]')

figure(2)
plot(t/60,Ys(:,2)-273.15,'r');
yline(Air_Temperature-273.15,'g')
title('Temperature of the liquid in time')
ylabel('Temperature of the liquid [C]')
title('time [min]')
legend('Liquid Temperature','Air Temperature')

function out = odefun(t,inputs,Water_MassDensity,MugDiameter,Steam_DynamicViscosity...
    ,Steam_SpecificHeat,Air_ThermalConductivity,Steam_Diffusivity...
    ,CeramicMug_ThermalConductivity,Water_HeatOfVaporisation,Water_SpecificHeat,Air_RelativeHumidity...
    ,Air_Temperature,Pressure,Water_MolarMass,Air_MolarMass,Mug_Thickness)

FluidMass = inputs(1);
Temperature = inputs(2);

FluidHeight = FluidMass/(Water_MassDensity*pi*MugDiameter^2/4);
VapourPressure = @(T) 10^(4.6543 - 1435.264/(T - 64.848 ))*10^5;
A_lat = pi*MugDiameter*FluidHeight;
A_top = pi*MugDiameter^2/4;
A_bottom = A_top;

GravityEarth = 9.81;
CharLength = MugDiameter;

Steam_MassDensityInterface = Water_MolarMass*VapourPressure(Temperature)/(8.314*Temperature);
PartialDensityH2O_interface = Steam_MassDensityInterface;
Steam_MaximumInAirMolarFraction = VapourPressure(Air_Temperature)/Pressure;
Steam_InAirMolarFraction = Steam_MaximumInAirMolarFraction*Air_RelativeHumidity;
HumidAir_AverageMolarMass = Steam_InAirMolarFraction*Water_MolarMass + (1-Steam_InAirMolarFraction)*Air_MolarMass;
Steam_InAirMassFraction = Steam_InAirMolarFraction*Water_MolarMass/HumidAir_AverageMolarMass;
HumidAir_BulkDensity = HumidAir_AverageMolarMass*Pressure/(8.314*Temperature);

PartialDensityH2O_bulk = HumidAir_BulkDensity*Steam_InAirMassFraction;


Grashoff = GravityEarth*CharLength^3/(Steam_DynamicViscosity/Steam_MassDensityInterface)^2*(HumidAir_BulkDensity/Steam_MassDensityInterface - 1);
Prandtl = Steam_DynamicViscosity*Steam_SpecificHeat/Air_ThermalConductivity;
Schmidt = Steam_DynamicViscosity/(Steam_MassDensityInterface*Steam_Diffusivity);
Rayleigh = Prandtl*Grashoff;

Sherwood = 0.0249*Grashoff^(2/5)*Schmidt^(7/15)*(1+0.494*Schmidt^(2/3))^(-2/5);
%Nusselt =  0.0249*Grashoff^(2/5)*Prandtl^(7/15)*(1+0.494*Prandtl^(2/3))^(-2/5);

MassTransferCoeff = Sherwood*Steam_Diffusivity/CharLength;
%HeatTransferCoeff = Sherwood*Steam_Diffusivity/CharLength;

HeatTransferCoeff_LateralExternal = (0.825 + 0.387*Rayleigh^(1/6)/((1+0.492/Prandtl^(9/16))^(8/27)))^2;


HeatTransferCoeff_LateralLayerCeramics = CeramicMug_ThermalConductivity/Mug_Thickness;
HeatTransferCoeff_Lateral = (1/HeatTransferCoeff_LateralExternal + 1/HeatTransferCoeff_LateralLayerCeramics)^-1;
HeatTransferCoeff_Bottom = CeramicMug_ThermalConductivity/(Mug_Thickness+3*Mug_Thickness);

ConductiveHeatExchanged = (A_lat*HeatTransferCoeff_Lateral + A_bottom*HeatTransferCoeff_Bottom)*(Temperature-Air_Temperature);


dmdt = -MassTransferCoeff*A_top*(PartialDensityH2O_interface - PartialDensityH2O_bulk);
%dmdt = 0;
dTdt = (- ConductiveHeatExchanged + dmdt*Water_HeatOfVaporisation)/(FluidMass*Water_SpecificHeat);
if FluidMass<0
    dmdt = 0;
else
end

out = [dmdt,dTdt]';
end



