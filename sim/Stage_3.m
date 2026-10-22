clear, clc, close all;
format shortE
%% Inputs - Stage 1
T               = 1939.3;           %Thrust in Newtons
m               = 6.369;            %mass in kg
g               = 9.81;             %grav accel in m/s^2
t_burn          = 2.232;            %burn time in seconds
dt              = 0.00025;            %change in time, gives exactly 2232 steps
n_burn          = round(t_burn / dt);   %steps in the burn
v_handcalc      = 657.728;          %hand-calculated velocity, m/s
h_handcalc      = 734.024;          %hand-calculated height, m
%% Inputs - Stage 2
d               = 0.079248;         %diameter of reference area S; in m
Cd              = 0.44;             %assumed constant drag coefficient
rho_1             = 1.225;            % density of air at sea level, kg/m^3
S               = (pi*(d^2))/4;      %reference area, in m^2
kd              = 0.5*rho_1*S*Cd;      %quadratic drag coefficient
D_upper         = kd*(v_handcalc^2); %Upperbound Drag force, uses no-drag burnout v; in N
t_max           = 100;               %generous on time to apogee; in seconds
N               = round(t_max / dt) +1; %Steps in the entire sim
a               = zeros(1,N);      % acceleration array, replaces stage 1 accel
%% Inputs - Stage 3
launchAlt       = 0;                 % launch site elevation above sea level, m
useAtm = true; 
%creating arrays
v = zeros(1, N);
h = zeros(1, N);
t = (0:N-1)*dt; %time array in seconds
D_hist = zeros(1, N);
rho_hist = zeros(1, N);

%% for loop using Euler Method
for i=(2:N)
    if i-1 <= n_burn
        Thr = T;
    else, Thr = 0;
    end
    z = launchAlt + h(i-1);
    if ~useAtm 
        rho = rho_1;
    else
        rho = stdatm(z);
    end
    D = 0.5.*rho.*S.*Cd.*(v(i-1)*abs(v(i-1)));
    D_hist(i-1) = D;
    rho_hist(i-1) = rho;
    a(i-1) = (Thr - m*g - D)/m;
    v(i) = v(i-1) + a(i-1) * dt;
    h(i) = h(i-1) + v(i-1) * dt;
    if v(i) < 0, break,
    end
end
%Trimming Values
v = v(1:i);
t = t(1:i);
h = h(1:i);
a = a(1:i-1);
D_hist = D_hist(1:i-1);
rho_hist = rho_hist(1:i-1);

[D_peak, k] = max(D_hist);
[h_apo, idx] = max(h);

fprintf('Burnout velocity is %.4f m/s at an altitude of %.4f m\n', v(n_burn+1), h(n_burn+1))
if i == N
    fprintf('Reached %.4f before apogee; increase t_max\n', t_max);
else
    fprintf(['Apogee: %.4f m | %.4f ft\n' ...
        'Time to Apogee: %.4f s\n' ...
        'Peak Drag: D = %.4f N | Upper Bound: %.4f N\n' ...
        'Time at Max Drag: t = %.4f s\n'], h_apo, h_apo/0.3048, t(idx),D_peak, D_upper, t(k))
end

%% creating figures
set(0, 'DefaultFigureVisible', 'on'); % TURNS OFF/ON figures while coding
%acceleration over time
figure
plot(t(1:i-1), a, 'g', Linewidth = 2);
ylabel('Acceleration (m/s^2)')
xlabel('Time (s)')
title('Acceleration over Time')
grid on
%velocity over time
figure
plot(t,v,'b',Linewidth=1)
ylabel ('Velocity (m/s)')
xlabel ('Time (s)')
title ('Velocity over Time')
grid on
%Height over Time
figure
plot(t,h,'r',Linewidth=1)
ylabel ('Height (m)')
xlabel ('Time (s)')
title ('Height over Time')
grid on
% Density over Altitude
figure
plot(rho_hist, h(1:i-1), 'c', Linewidth = 1)
ylabel('Height (m)')
xlabel('Density (kg/m^3)')
title('Density over Altitude')
grid on
% Drag over Time
figure
plot(t(1:i-1), D_hist, 'm', Linewidth = 1)
xlabel('Time (s)')
ylabel('Drag (N)')
title('Drag over Time')
grid on
