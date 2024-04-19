%% This function is to simulate the system of a single-phase inverter 
%% operating under sine-triangle PWM with the load is an rl load with 
%% r = 0.5 ohm and l = 1 mH and the desired ac voltage has a fundamental
%% frequency of 400Hz and a phase angle phi = -pi/2. We assume DC voltage
%% is 100V and switching frequency is 7600Hz.
clear
%% Initial Parameters
r = 1; % resistive load
l = 1e-3; % inductive load
phi = (-pi / 2) * 180 / pi; % phase angle
f = 400; % fundamental frequency
fsw = 7600; % switching frequency
T = 1 / f; % fundamental period
Vdc = 100; % DC voltage
w_ac = 2 * pi * f; % omega
w_sw = 2 * pi * fsw; % omega switch
dt = 1e-7; % time step
tend = 2*T; % end of loop
t = 0:dt:tend; % time array
t_len = length(t);

tri = 1 / 2; % tri function intial value
for i = 1:30
    tri = tri + (2 * cos(pi * i) - cos(2 * pi * i) - 1) / ((pi ^ 2) ...
            * (i ^ 2)) * cos(w_sw * i * t);
end
tri = 2 * tri - 1;
m = 1; 

%% Components Arrays of Values
Vac = zeros(1, t_len);
Idc = zeros(1, t_len);
Iac = zeros(1, t_len);
Id1 = zeros(1, t_len);
Id2 = zeros(1, t_len);
Id3 = zeros(1, t_len);
Id4 = zeros(1, t_len);
It1 = zeros(1, t_len);
It2 = zeros(1, t_len);
It3 = zeros(1, t_len);
It4 = zeros(1, t_len);
theta_ac = zeros(1, t_len);
Vd1 = zeros(1, t_len);
Vd2 = zeros(1, t_len);
Vd3 = zeros(1, t_len);
Vd4 = zeros(1, t_len);
Vt1 = zeros(1, t_len);
Vt2 = zeros(1, t_len);
Vt3 = zeros(1, t_len);
Vt4 = zeros(1, t_len);
d = zeros(1, t_len); % duty cycle array

%% while loop for simulation
k = 1; % while loop counter
while (t(k) < tend)
    theta_ac(k + 1) = (w_ac * t(k)) * 180 / pi;
    theta_ac(k + 1) = mod(theta_ac(k + 1), 360);
    d(k) = m * cosd(theta_ac(k) + phi);
    %% Switching logic
    if (d(k) >= tri(k))
        T2 = 0;
        T3 = 0;
        T1 = 1;
        T4 = 1;
    elseif (d(k) < tri(k))
        T2 = 1;
        T3 = 1;
        T1 = 0;
        T4 = 0;
    end
    %% AC Voltage and Current
    Vac(k + 1) = (2 * T1 * T4 - 1) * Vdc;
    Iac(k + 1) = Iac(k) + dt * ((Vac(k) - r * Iac(k)) / l);
    %% Diodes and Transistors Voltages
    Vt1(k) = (1 - T1) * Vdc;
    Vt4(k) = Vt1(k);
    Vt2(k) = (1 - T2) * Vdc;
    Vt3(k) = Vt2(k);
    Vd1(k) = -T1 * Vdc;
    Vd4(k) = Vd1(k);
    Vd2(k) = -T2 * Vdc;
    Vd3(k) = Vd2(k);
    %% Diodes and Transistors Currents
    if (T1 == 1 && T4 == 1)
        if (Iac(k) >= 0) 
            It1(k) = Iac(k);
            It4(k) = It1(k);
            Id1(k) = 0;
            Id4(k) = 0;
        else
            It1(k) = 0;
            It4(k) = 0;
            Id1(k) = -Iac(k);
            Id4(k) = Id1(k);
        end
    elseif (T2 == 1 && T3 == 1)
        if (Iac(k) >= 0) 
            Id2(k) = Iac(k);
            Id3(k) = Iac(k);
            It2(k) = 0;
            It3(k) = 0;
        else
            Id2(k) = 0;
            Id3(k) = 0;
            It2(k) = -Iac(k);
            It3(k) = -Iac(k);
        end
    end
    Idc(k) = It1(k) + It2(k) - Id1(k) - Id2(k);
    % Time array and steps update
    t(k + 1) = t(k) + dt;
    k = k + 1;
end

%% Calculate frequency spectrum
f_spec = f:f:fsw;
N = length(f_spec);
[avg, ak, bk, rcon, err] = fourseries(t, Vac, T, N);
disp(err)
disp(avg)

%% Plot
% Note: the indexing seen here in every arrays is to compensate for the mod
% function used earlier to calculate and update theta_ac. This function,
% after hitting 360, returns straight back to zero, which results in an
% overlap of many periods, which affects the integrity of the analysis.
% Index 25003 is where theta_ac returns back to zero.
figure(1)
subplot(2,1,1)
plot(theta_ac(25003:end), Vac(25003:end)) 
ylim([-120 120])
grid on
title('Vac(V) vs \theta_{AC}(degree)')
subplot(2,1,2)
plot(theta_ac(25003:end), Iac(25003:end))
ylim([-60 60])
grid on
title('Iac(A) vs \theta_{AC}(degree)')

figure(2)
subplot(2,2,1)
plot(theta_ac(25003:end), Vd1(25003:end))
ylim([-110 10])
grid on
title('Vd1(V) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(25003:end), Vd2(25003:end))
ylim([-110 10])
grid on
title('Vd2(V) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(25003:end), Vd3(25003:end))
ylim([-110 10])
grid on
title('Vd3(V) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(25003:end), Vd4(25003:end))
ylim([-110 10])
grid on
title('Vd4(V) vs \theta_{AC}(degree)')

figure(3)
subplot(2,2,1)
plot(theta_ac(25003:end), Vt1(25003:end))
ylim([-10 110])
grid on
title('Vt1(V) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(25003:end), Vt2(25003:end))
ylim([-10 110])
grid on
title('Vt2(V) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(25003:end), Vt3(25003:end))
ylim([-10 110])
grid on
title('Vt3(V) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(25003:end), Vt4(25003:end))
ylim([-10 110])
grid on
title('Vt4(V) vs \theta_{AC}(degree)')

figure(4)
subplot(2,2,1)
plot(theta_ac(25003:end), Id1(25003:end))
grid on
title('Id1(A) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(25003:end), Id2(25003:end))
grid on
title('Id2(A) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(25003:end), Id3(25003:end))
grid on
title('Id3(A) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(25003:end), Id4(25003:end))
grid on
title('Id4(A) vs \theta_{AC}(degree)')

figure(5)
subplot(2,2,1)
plot(theta_ac(25003:end), It1(25003:end))
grid on
title('It1(A) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(25003:end), It2(25003:end))
grid on
title('It2(A) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(25003:end), It3(25003:end))
grid on
title('It3(A) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(25003:end), It4(25003:end))
grid on
title('It4(A) vs \theta_{AC}(degree)')

figure(6)
plot(theta_ac(25003:end), Idc(25003:end))
grid on
title('Idc(A) vs \theta_{AC}(degree)')

figure(7)
plot(f_spec, ak, 'o')
hold on
plot(f_spec, bk, '*')
title('Frequency Spectrum')
ylim([-20 120])
legend('ak','bk')
grid on
hold off

figure(8)
plot(t, d, t, tri)
legend('Duty Cycle','Triangle Wave')
title('Duty Cycle and Triangle Wave vs Time')
grid on