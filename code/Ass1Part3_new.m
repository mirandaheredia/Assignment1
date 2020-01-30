%% Part 3 - Enhancements
%% 
%Question 3
%The purpose of this code is to model the electrons as particles with the
%effective mass using Monte-Carlo Modelling
%The scattering is in place here but the Mean free path and temp changes
%due to scatters
%This code contains the electron denisity map
%Miranda Heredia 100996160

clear
close all

%rng(101);

%Constants given
width = 200e-9;
height = 100e-9;
iterations = 10;
Noelectron = 100;
K = 1.38064e-23;
T = 300;
mo = 9.1093856e-31;
m = 0.26*mo;
vth = sqrt(2*K*T/m);
Time_step = .002e-12;

%Initialize the velocities and positions with an empty array
Vel_x = zeros(Noelectron,1);
Vel_y = zeros(Noelectron,1);
Pos_x = zeros(Noelectron,1);
Pos_y = zeros(Noelectron,1);


%randomize initial positions
Pos_x(:,1) = width*rand(Noelectron,1);
Pos_y(:,1) = height*rand(Noelectron,1);

%Initialize boundaries for Boxes
xBox = Pos_x >80e-9 & Pos_x<120e-9; %note same box boundary will exist for x region
UpperBox = xBox & Pos_y > 60e-9;
LowerBox = xBox & Pos_y< 40e-9;
Inside = UpperBox | LowerBox;

WidthBox = 40e-9;



%%
%if the particle starts inside, move outside of the box
%movement is not random and should be
for j=1:(length(Inside))
    
    if(Inside(j))
        Pos_x(j) = Pos_x(j) - WidthBox;
    end
end


%Plots the initial positions of the particles
figure(1)
plot(Pos_x,Pos_y, '.')
% axis ([0 200e-9 0 100e-9])
hold on
rectangle('Position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7])
rectangle('Position',[0.8e-7 0 0.4e-7 0.4e-7])
title('Initial Positions');
hold off

%Creates random distribution for velocites
for n=1:Noelectron
    Vel_x(n,1)= randn(1,1)*vth;
    Vel_y(n,1) = randn(1,1)*vth;
end


%Parameters set to plot temperature average
timesteps = 100;
T_avg = zeros(timesteps,1);
time = 1:timesteps;



%PScattering
colorarray = rand(Noelectron,1);
p_scatter = 1 - exp(-Time_step/0.2e-12);
for n= 1:timesteps
    Pos_x_old = Pos_x;
    Pos_y_old = Pos_y;

    random = rand(Noelectron,1);
    
    %all electrons with higher probabilities
    new = random < p_scatter;
    
    %all electrons with lower probabilities
    new2 = random >= p_scatter;
    

    rand_v_x = zeros(Noelectron,1);
    rand_v_y = zeros(Noelectron,1);
    
    for i = 1:1:Noelectron
     r1 = randi([1 Noelectron], 1,1);
     r2 = randi([1 Noelectron], 1,1);
        rand_v_x(i,1) = Vel_x(r1,1);
        rand_v_y(i,1) = Vel_y(r2,1);
    end
    %all electrons with lower probabilities will stay the same
       Vel_x = Vel_x.*new2;
       Vel_y = Vel_y.*new2;

       rand_v_x=rand_v_x.*new;
       rand_v_y=rand_v_y.*new;

       Vel_x = Vel_x+rand_v_x;
       Vel_y = Vel_y+rand_v_y;

       Pos_x = Pos_x + Vel_x*Time_step;
       Pos_y = Pos_y + Vel_y*Time_step;
    
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
    
    %Rectangle Boundary Conditions (bottleneck)
    idx = (Pos_x < 1.2e-7) & (Pos_x> 0.8e-7);
    idy = Pos_y >0.6e-7 | Pos_y < 0.4e-7;
    
    inBox = idx & idy;
    
    OutX = Pos_x_old < 0.8e-7 | Pos_x_old >1.2e-7; %not used
    BetweenY = Pos_y_old > 0.4e-7 & Pos_y_old < 0.6e-7;
    
    Vel_y(inBox & BetweenY) = -1*Vel_y(inBox & BetweenY);
    Vel_x(inBox & ~BetweenY) = -1*Vel_x(inBox & ~BetweenY);
    
    

    v_avg = mean((Vel_x.^2)+(Vel_y.^2));
    T_avg(n) =(1/2)* (m*(v_avg))*(1/K);
    
    
    % mean free path 
     mfp = (10^-15)*(v_avg);
    
    %%
    %2-D plot of particle trajectories
    figure(3)
    scatter(Pos_x,Pos_y,5,colorarray,'o');
    axis([0 200*10^-9 0 100*10^-9])
    title(['The mean free path is ', num2str(mfp)]);
    hold on
    rectangle('Position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7])
    rectangle('Position',[0.8e-7 0 0.4e-7 0.4e-7])
    
    
    %%
    % Average Temperature Plot
    figure(4)
    plot(time,T_avg)
    title(['The Average Temperature is ', num2str(T_avg(n))]);
    
 
    
    %average thermal velocity
   
end
%%
%Electron Density Plot
figure(5)
hist3([Pos_x Pos_y],'CdataMode','auto');
title('Electron Density Map of Final Positions')
colorbar
view(2)


