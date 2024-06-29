clear
clc
close all




ICs = [0,300/3.6];

%[BestCombination,MaxSpeed] = ga(@(B,Time_engamentElectricMotor)optimiserFCN(B,Time_engamentElectricMotor),2,[],[],[],[],[0.001,0],[1000,120]);

%BestCombination = ga(@(B,Time_engamentElectricMotor)optimiserFCN(B,Time_engamentElectricMotor),2,[],[],[],[],[0.001,0.001],[1000,120],[]);

BestCombination = ga(@optimiserFCN,2,[],[],[],[],[0.1,0.1],[1000,100]);

B = BestCombination(1)

Time_engamentElectricMotor = BestCombination(2)

%Disp(MaxSpeed)

[timeCoor,Ys_solution] = ode23s(@(t,Ys)odefun(t,Ys,B,Time_engamentElectricMotor),[0,150],ICs);


figure(1)

plot(timeCoor,Ys_solution(:,2)*3.6)
title('Speed-time relationship')
xlabel('time [s]')
ylabel('speed [km/h]')


figure(2)
plot(Ys_solution(:,1),Ys_solution(:,2))
title('Speed-space relationship')
xlabel('spatial coordinate [m]')
ylabel('speed [km/h]')



function out = optimiserFCN(InputVars)

ICs = [0,300/3.6];
B = InputVars(1);
Time_engamentElectricMotor = InputVars(2);

[timeCoor,Ys_solution] = ode23s(@(t,Ys)odefun(t,Ys,B,Time_engamentElectricMotor),[0,300],ICs);

 [MaxSpeed,Index] = max(Ys_solution(:,2)*3.6);

out = 1/abs(MaxSpeed);

if Ys_solution(Index,1) > 8000
    out = out + 10;
else
end

end


function out = odefun(t,Ys,B,Time_engamentElectricMotor)

SpatialCoor = Ys(1);

Velocity = Ys(2);

DragCoefficient = 0.333; % -, reddit source

Mass = 1990; % kg

PowerICE = 735499; % W

FrontalArea = 2.142; %m2

PowerElectricMotor_Max_FullCharge = 800*735499/1000; % W

PowerElectricMotor_Max = PowerElectricMotor_Max_FullCharge*(1-0.01*(t-Time_engamentElectricMotor));

if PowerElectricMotor_Max < 0
    PowerElectricMotor_Max = 0;
elseif t < Time_engamentElectricMotor
    PowerElectricMotor_Max = PowerElectricMotor_Max_FullCharge;
end

PowerElectricMotor = PowerElectricMotor_Max*(B*(t-Time_engamentElectricMotor));

if PowerElectricMotor > PowerElectricMotor_Max
    PowerElectricMotor = PowerElectricMotor_Max;
elseif PowerElectricMotor < 0
    PowerElectricMotor = 0;
end
    

DensityAir = 1.2250; %kg/m3

DriveTrainEfficiency = 0.85; %  -

TotalPower = PowerElectricMotor + PowerICE;

Velocity_kmh = Velocity/3.6;

RollingFrictionCoefficient = 0.005 + (1 / 3.44738)*(0.01 + 0.0095*(Velocity_kmh / 100)^2); %engineering toolbox

RollingFriction = RollingFrictionCoefficient*Mass*9.81;

DragForce = 1/2*FrontalArea*DragCoefficient*DensityAir*Velocity^2;

dSpaceCoordt = Velocity;

dVelocitydt = TotalPower*DriveTrainEfficiency/(Mass*Velocity) - DragForce/Mass - RollingFriction/Mass;


out = [dSpaceCoordt dVelocitydt]';
end