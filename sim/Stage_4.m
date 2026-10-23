% Stage 4: thrust from a RASP curve file,
% mass drops as propellant burns
clear, clc, close all;
format shortE
%% Inputs
m_wet           = 6.369;            %mass in kg
g               = 9.81;             %grav accel in m/s^2
dt              = 0.00025;            %change in time
d               = 0.079248;         %diameter of reference area S; in m
Cd              = 0.44;             %assumed constant drag coefficient
rho_1           = 1.225;            % density of air at sea level, kg/m^3
S               = (pi*(d^2))/4;      %reference area, in m^2
t_max           = 100;               %generous on time to apogee; in seconds
N               = round(t_max / dt) +1; %Steps in the entire sim
a               = zeros(1,N);      % acceleration array
launchAlt       = 0;                 % launch site elevation above sea level, m
useAtm          = true;
%% Stage 4 Inputs
motorFile       = 'data/AeroTech_L1940X.eng';
useMassDrop     = true;
%% creating arrays
v = zeros(1, N);
h = zeros(1, N);
t = (0:N-1)*dt; %time array in seconds
D_hist = zeros(1, N);
rho_hist = zeros(1, N);
[tc, Fc, m_prop] = readEng(motorFile);
t_bo = tc(end);
F_grid = interp1(tc, Fc, t, "linear", 0); % gives thrust at every time, and 0 after burnout
%% Precomputing Mass
I_cum = cumtrapz(tc, Fc);
frac_c = I_cum / I_cum(end);
frac_grid = interp1(tc, frac_c, t, 'linear', 1);
if useMassDrop
    m_grid = m_wet - m_prop * frac_grid;
else
    m_grid = m_wet * ones(size(t));
end
%% for loop using Euler Method
for i=(2:N)
    Thr = F_grid(i-1);
    z = launchAlt + h(i-1);
    if ~useAtm
        rho = rho_1;
    else
        rho = stdatm(z);
    end
    D = 0.5.*rho.*S.*Cd.*(v(i-1)*abs(v(i-1)));
    D_hist(i-1) = D;
    rho_hist(i-1) = rho;
    a(i-1) = (Thr - (m_grid(i-1))*g - D)/(m_grid(i-1));
    if h(i-1) <=0 && a(i-1)<0
        a(i-1) = 0;
    end
    v(i) = v(i-1) + a(i-1) * dt;
    h(i) = h(i-1) + v(i-1) * dt;
    if v(i) < 0, break,
    end
end
%% Trimming Values
v = v(1:i);
t = t(1:i);
h = h(1:i);
a = a(1:i-1);
D_hist = D_hist(1:i-1);
rho_hist = rho_hist(1:i-1);
F_grid = F_grid(1:i);
m_grid = m_grid(1:i);

[D_peak, k] = max(D_hist);
[h_apo, idx] = max(h);
i_bo = find(t >= t_bo, 1);
[v_max, i_vmax] = max(v);
%% Prints
fprintf(['Burnout velocity is %.4f m/s at an altitude of %.4f m, happening at %.4fs\n ' ...
    'Max Velocity: %.4f | Max Acceleration: %.4f'], v(i_bo), h(i_bo), t_bo, max(v), max(a))
if i == N
    fprintf('Reached %.4f before apogee; increase t_max\n', t_max);
else
    fprintf(['Apogee: %.4f m | %.4f ft\n' ...
        'Time to Apogee: %.4f s\n' ...
        'Peak Drag: D = %.4f N\n' ...
        'Time at Max Drag: t = %.4f s\n'...
        'Propellant Mass: %.4f kg\n', ...
        'Burnout Mass: %.4f\n'], ...
        h_apo, h_apo/0.3048, t(idx),D_peak, t(k), m_prop, m_grid(i_bo))


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
% Thrust over time
figure
plot(t, F_grid, 'r', Linewidth = 1)
xlabel('Time (s)')
ylabel('Thrust (N)')
title('Thrust over Time')
grid on
% Mass over time
figure
plot(t, m_grid, 'w', Linewidth = 1)
xlabel('Time (s)')
ylabel('Mass (kg)')
title('Mass over Time')
grid on