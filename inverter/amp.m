%% amp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This function is to simulate the system of a single-phase inverter 
operating under sinetriangle PWM. this function varies m values from 0, 
0.2, 0.4, 1.0, 1.3, 1.6, 2.0, 3.0, 4, 5 then plot amplitude as a function
of duty cycle
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Initial Parameters
r = 1; % resistive load
l = 1e-3; % inductive load
phi = (-pi / 2) * 180 / pi; % phase angle
f = 400; % fundamental frequency
T = 1 / f; % fundamental period
w_ac = 2 * pi * f; % omega
fsw = 7600; % switching frequency
Tsw = 1 / fsw;
w_sw = 2 * pi * fsw; % omega switch

Vdc = 100; % DC voltage
dt = 1e-7; % time step
tend = 2*T; % end of loop
t = 0:dt:tend; % time array
t_len = length(t);
m_array = [0, 0.2, 0.4, 1.0, 1.3, 1.6, 2.0, 3.0, 4, 5];
m_len = length(m_array);

tri = 1 / 2; % tri function intial value
for i = 1:30
    tri = tri + (2 * cos(pi * i) - cos(2 * pi * i) - 1) / ((pi ^ 2) ...
            * (i ^ 2)) * cos(w_sw * i * t);
end
tri = 2 * tri - 1;
m = 1; 

%% Components Arrays of Values
Vac = zeros(1, t_len);
theta_ac = zeros(1, t_len);
d = zeros(1, t_len);
amplitude = zeros(1, m_len);

%% Procedure
for i = 1:m_len
    %% while loop for simulation
    k = 1; % while loop counter
    while (t(k) < tend)
        theta_ac(k + 1) = (w_ac * t(k)) * 180 / pi;
        theta_ac(k + 1) = mod(theta_ac(k + 1), 360);
        d(k) = m_array(i) * cosd(theta_ac(k) + phi);
        %% Switching logic
        if (d(k) >= tri(k))
            T1 = 1;
            T4 = 1;
        else
            T1 = 0;
            T4 = 0;
        end
        %% AC Voltage 
        Vac(k + 1) = (2 * T1 * T4 - 1) * Vdc;
        % Time array and steps update
        t(k + 1) = t(k) + dt;
        k = k + 1;
    end
    [avg,ak,bk,rcon,err] = fourseries(t, Vac, T, 10);
    amplitude(i) = sqrt((ak(1)) ^ 2 + (bk(1)) ^ 2);
end

%% Plot
figure(1)
plot(m_array, amplitude);
title('AC Voltage Amplitudes vs Duty Cycle')
xlabel('m values')
ylabel('AC Voltage Amplitudes')
grid on