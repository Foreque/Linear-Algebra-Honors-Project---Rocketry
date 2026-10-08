clear,clc, close all
format short
%Inputs
T = 1939.3; %Thrust in Newtons
m = 6.369; %mass in kg
g = 9.81; %grav accel in m/s^2
BT = 2.232; %burn time in seconds
dt = 0.001; %change in time, gives exactly 2232 steps
N = round(BT / dt);%steps in the sim
t_placeholder = (0:N); %creating a placeholder so array t can be created
a = (T - m*g)/m; %solving acceleration; derived from F=ma at t=0+
v_handcalc = 657.728; %hand-calculated velocity, m/s
h_handcalc = 734.024; %hand-calculated height, m
%creating arrays
v = zeros(1, N+1);
h = zeros(1, N+1);
t = t_placeholder.*dt; %time array in seconds

%for loop with Eulers Method added
for k=(2:N+1)
    v(k) = v(k-1) + a * dt;
    h(k) = h(k-1) + v(k-1) * dt;
end

%% creating figures
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


%% Print out final values

fprintf(['Final velocity is %.4f m/s. My paper calculations gave %.4f m/s, which MATLAB agrees to the ' ...
'thousandths place.\n'], v(end), v_handcalc) 
fprintf(['Final Height is %.4f m. My paper calculations gave %.4f m, a difference of %.4f.\n' ...
'I believe this is due to Eulers Method increasing velocity with each step, and my hand calculation\n' ...
'is using the true values; no approximation using Eulers method.\n'], h(end), h_handcalc, h_handcalc - h(end))