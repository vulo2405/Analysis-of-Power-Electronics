%% switch180.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This function is to simulate the system of a single-phase inverter 
operating under 180 switching with the load is an rl load with 
r = 0.5 ohm and l = 1 mH and the desired fundamental component
of AC voltage is Vas = -200/pi*cos(theta_ac) where theta_ac = 120*pi*t
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Initial Parameters
r = 0.5; % resistive load
l = 1e-3; % inductive load
f = 60; % fundamental frequency
T = 1 / f; % fundamental period
Vdc = 50; % Vdc
w = 2 * pi * f;
dt = 1e-7; % time step
tend = 2*T; % end of loop
t = 0:dt:tend; % time array
t_len = length(t);

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

%% while loop for simulation
k = 1; % while loop counter
while (t(k) < tend)
    theta_ac(k + 1) = w * t(k) * 180 / pi;
    theta_ac(k + 1) = mod(theta_ac(k + 1), 360);
    %% Switching logic
    if (theta_ac(k) >= 0 && theta_ac(k) < 90) || (theta_ac(k) > 270 && theta_ac(k) <= 360)
        T2 = 1;
        T3 = 1;
        T1 = 0;
        T4 = 0;
    elseif (theta_ac(k) >= 90 && theta_ac(k) <= 270)
        T2 = 0;
        T3 = 0;
        T1 = 1;
        T4 = 1;
    end
    %% AC Voltage and Current
    Vac(k + 1) = (2 * T1 * T4 - 1) * Vdc;
    Iac(k + 1) = Iac(k) + dt * (Vac(k) - r * Iac(k)) / l;
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
            Id4(k) = -Iac(k);
        end
    else
        if (Iac(k) >= 0) 
            It2(k) = 0;
            It3(k) = 0;
            Id2(k) = Iac(k);
            Id3(k) = Iac(k);
        else
            It2(k) = -Iac(k);
            It3(k) = -Iac(k);
            Id2(k) = 0;
            Id3(k) = 0;
        end
    end
    Idc(k) = It1(k) + It2(k) - Id1(k) - Id2(k);
    % Time array and steps update
    t(k + 1) = t(k) + dt;
    k = k + 1;
end

%% Calculate frequency spectrum
f_spec = f:f:1200;
N = length(f_spec);
[avg, ak, bk, rcon, err] = fourseries(t, Vac, T, N);
disp(err)

%% Plot
% Note: the indexing seen here in every arrays is to compensate for the mod
% function used earlier to calculate and update theta_ac. This function,
% after hitting 360, returns straight back to zero, which results in an
% overlap of many periods, which affects the integrity of the analysis.
% Index 166669 is where theta_ac returns back to zero.
figure(1)
subplot(2,1,1)
plot(theta_ac(166669:end), Vac(166669:end)) 
ylim([-80 80])
grid on
title('Vac(V) vs \theta_{AC}(degree)')
subplot(2,1,2)
plot(theta_ac(166669:end), Iac(166669:end))
ylim([-120 120])
grid on
title('Iac(A) vs \theta_{AC}(degree)')

figure(2)
subplot(2,2,1)
plot(theta_ac(166669:end), Vd1(166668:end))
ylim([-80 20])
grid on
title('Vd1(V) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(166669:end), Vd2(166668:end))
ylim([-80 20])
grid on
title('Vd2(V) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(166669:end), Vd3(166668:end))
ylim([-80 20])
grid on
title('Vd3(V) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(166669:end), Vd4(166668:end))
ylim([-80 20])
grid on
title('Vd4(V) vs \theta_{AC}(degree)')

figure(3)
subplot(2,2,1)
plot(theta_ac(166669:end), Vt1(166668:end))
ylim([-20 80])
grid on
title('Vt1(V) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(166669:end), Vt2(166668:end))
ylim([-20 80])
grid on
title('Vt2(V) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(166669:end), Vt3(166668:end))
ylim([-20 80])
grid on
title('Vt3(V) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(166669:end), Vt4(166668:end))
ylim([-20 80])
grid on
title('Vt4(V) vs \theta_{AC}(degree)')

figure(4)
subplot(2,2,1)
plot(theta_ac(166669:end), Id1(166668:end))
grid on
title('Id1(A) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(166669:end), Id2(166668:end))
grid on
title('Id2(A) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(166669:end), Id3(166668:end))
grid on
title('Id3(A) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(166669:end), Id4(166668:end))
grid on
title('Id4(A) vs \theta_{AC}(degree)')

figure(5)
subplot(2,2,1)
plot(theta_ac(166669:end), It1(166668:end))
grid on
title('It1(A) vs \theta_{AC}(degree)')
subplot(2,2,2)
plot(theta_ac(166669:end), It2(166668:end))
grid on
title('It2(A) vs \theta_{AC}(degree)')
subplot(2,2,3)
plot(theta_ac(166669:end), It3(166668:end))
grid on
title('It3(A) vs \theta_{AC}(degree)')
subplot(2,2,4)
plot(theta_ac(166669:end), It4(166668:end))
grid on
title('It4(A) vs \theta_{AC}(degree)')

figure(6)
plot(theta_ac(166669:end), Idc(166668:end))
grid on
title('Idc(A) vs \theta_{AC}(degree)')

figure(7)
plot(f_spec, ak, 'o')
hold on
plot(f_spec, bk, '+')
title('Frequency Spectrum')
legend('ak','bk')
grid on
hold off