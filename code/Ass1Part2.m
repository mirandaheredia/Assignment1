%% Part 2 Collisions with Mean Free Path


%Constants given
length = 200e-9;
height = 100e-9;
iterations = 10;
Noelectron = 100;
K = 1.38064e-23;
T = 300;
mo = 9.1093856e-31;
m = 0.26*mo;
vth = sqrt(2*K*T/m);
Time_step = .002e-12;


Vel_x = zeros(Noelectron,1);
Vel_y = zeros(Noelectron,1);

Pos_x = zeros(Noelectron,1);
Pos_y = zeros(Noelectron,1);

Pos_x(:,1) = length*rand(Noelectron,1);
Pos_y(:,1) = height*rand(Noelectron,1);

for n=1:Noelectron
    Vel_x(n,1)= randn(1,1)*vth;
    Vel_y(n,1) = randn(1,1)*vth;
end
AvgV = sqrt(Vel_x.^2 + Vel_y.^2);
figure(1)
hist(Vel_x,Noelectron)
title('X component of distribution')
figure(2)
hist(Vel_y, Noelectron)
title('Y component of distribution')
figure(3)
hist(AvgV,Noelectron)

%initialize temperature vector
timesteps = 100;
T_avg_V = zeros(timesteps,1);


%PScattering
colorarray = rand(Noelectron,1);
p_scatter = 1 - exp(-Time_step/0.2e-12);
time = 1:timesteps;
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
   Pos_y = Pos_y +Vel_y*Time_step;
    
   %checking for boundary positions
    idLong = Pos_x>=length;
    Pos_x(idLong) = Pos_x(idLong) - length;
    Pos_x_old(idLong) = 0;
    
    idShort = Pos_x<=0;
    Pos_x(idShort) = Pos_x(idShort) + length;
    Pos_x_old(idShort) = length;
    
    %Check for y boundary and correct
       
    Vel_y(Pos_y>=height) = -1*Vel_y(Pos_y>=height);
    Vel_y(Pos_y<=0) = -1*Vel_y(Pos_y<=0);
    
    Pos_y(Pos_y>height) = height - (Pos_y(Pos_y>height)-height);
   
    %average thermal velocity
    v_avg = mean((Vel_x.^2)+(Vel_y.^2));
    T_avg_V(n) = (1/2)*(m*(v_avg))/K;
    
    
    %%
    % mean free path is 0.0127 m
    
    mfp = (0.2e-15)*(v_avg);
    %%
    %mean time between collisions is 1.2639e-05
    meantime = (Noelectron./sum(v_avg))*mfp;
    
    figure(4)
    scatter(Pos_x,Pos_y,3,colorarray);
    axis([0 200e-9 0 100e-9])
    title(['The mean free path is ', num2str(mfp)]);
    hold on
    
    figure(5)
    plot(time,T_avg_V)
    title(['The Average Temperature is ', num2str(T_avg_V(n))]);
   
end

