%Miranda Heredia 100996160
%% Part 1 - Electron Modelling

clear
close all
%Constants given
width = 200e-9;
height = 100e-9;
iterations = 100;
Noelectron = 100;
K = 1.38064e-23;
T = 300;
mo = 9.1093856e-31;
m = 0.26*mo;

%%
% The thermal velocity is 187 km/s calculated with 2 degrees of freedom
vth = sqrt(2*K*T/m);

%%
% The mean free path is 37.4 nm
mfp= vth*0.2e-12;

%array to hold info
Pos_x = zeros(Noelectron,1);
Pos_y = zeros(Noelectron,1);
Vel_x = zeros(Noelectron,1);
Vel_y = zeros(Noelectron,1);

%Initial array population

Pos_x(:,1) = width*rand(Noelectron,1);
Pos_y(:,1) = height*rand(Noelectron,1);


%Plot initial positions
figure(1)
plot(Pos_x,Pos_y, 'o')
title('Initial Positions')

%Arbitrary time step
Time_step = .02e-12;

%randomize initial velocities
angle = 2*pi*rand(Noelectron,1);
Vel_x = vth*cos(angle);
Vel_y = vth*sin(angle);

colorarray = rand(Noelectron,1);
timesteps = iterations;
T_avg_V = zeros(timesteps,1);
time = 1:timesteps;

for i = 1:iterations
    
    %mainitain old positions
    Pos_x_old = Pos_x;
    Pos_y_old = Pos_y;
    
    %Update the Positions
    Pos_x = Pos_x + Vel_x*Time_step;
    Pos_y = Pos_y +Vel_y*Time_step;
    
    %checking for boundary positions
    idLong = Pos_x>=width;
    Pos_x(idLong) = Pos_x(idLong) - width;
    Pos_x_old(idLong) = 0;
    
    idShort = Pos_x<=0;
    Pos_x(idShort) = Pos_x(idShort) + width;
    Pos_x_old(idShort) = width;
    
    %Check for y boundary and correct
       
    Vel_y(Pos_y>=height) = -1*Vel_y(Pos_y>=height);
    Vel_y(Pos_y<=0) = -1*Vel_y(Pos_y<=0);
    
    Pos_y(Pos_y>height) = height - (Pos_y(Pos_y>height)-height);
    %Pos_y(Pos_y<0) = height - (Pos_y(Pos_y<0)-height);
        
 
    figure(2)
    plot([Pos_x_old';Pos_x'], [Pos_y_old';Pos_y'],3,colorarray); 
    title('Electron Modelling')
   
    axis ([0 200e-9 0 100e-9])
    
    hold on
    pause(0.05)
    
    %Average Temperature Plot
    figure(3)
    Vel_avg_2 =mean((Vel_x.^2)+(Vel_y.^2));
    T_avg_V(i) = (1/2)*(1/K)*(m*(Vel_avg_2));
    figure(3)
    plot(time,T_avg_V);
    title(['The Average Temperature is ', num2str(T_avg_V(i))])
    
end
 


