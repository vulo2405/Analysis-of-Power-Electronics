%% sinetri3phase.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the system of a 3-phase inverter 
operating under sine-triangle PWM with the load is an rl load with 
r = 0.5 ohm and l = 1 mH and the desired ac current has a fundamental
component at 50Hz and phase a current is Ias = 20cos(theta_ac). 
We assume switching frequency is 3100Hz.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Initial Parameters
r = 1; % resistive load
l = 1e-3; % inductive load
phi_dis = 120; % displayed phase for each leg

f = 50; % fundamental frequency
T = 1 / f; % fundamental period
w_ac = 2 * pi * f; % omega
fsw = 3100; % switching frequency
w_sw = 2 * pi * fsw; % omega switch
T_sw = 1 / fsw; % switching period

dt = 1e-7; % time step
tend = 2*T; % end of loop
t = 0:dt:tend; % time array
t_len = length(t);

%% Calculate Vdc and phase shift angle
Vdc = 2 * sqrt(20 ^ 2 + (2 * pi) ^ 2); % Vdc derived from given Ias
phi = atand((-2 * pi) / 20); % Phase shift angle

%% Triangle function
tri = 1 / 2; % tri function intial value
for i = 1:100
    tri = tri + (2 * cos(pi * i) - cos(2 * pi * i) - 1) / ((pi ^ 2) ...
            * (i ^ 2)) * cos(w_sw * i * t);
end
tri = 2 * tri - 1;
m = 1; 

%% Components Arrays of Values
Vas = zeros(1, t_len);
Vbs = zeros(1, t_len);
Vcs = zeros(1, t_len);
Vag = zeros(1, t_len);
Vbg = zeros(1, t_len);
Vcg = zeros(1, t_len);

Idc = zeros(1, t_len);

Ias = zeros(1, t_len);
Ibs = zeros(1, t_len);
Ics = zeros(1, t_len);

Id1 = zeros(1, t_len);
Id2 = zeros(1, t_len);
Id3 = zeros(1, t_len);
Id4 = zeros(1, t_len);
Id5 = zeros(1, t_len);
Id6 = zeros(1, t_len);

It1 = zeros(1, t_len);
It2 = zeros(1, t_len);
It3 = zeros(1, t_len);
It4 = zeros(1, t_len);
It5 = zeros(1, t_len);
It6 = zeros(1, t_len);

theta_ac = zeros(1, t_len);

Vd1 = zeros(1, t_len);
Vd2 = zeros(1, t_len);
Vd3 = zeros(1, t_len);
Vd4 = zeros(1, t_len);
Vd5 = zeros(1, t_len);
Vd6 = zeros(1, t_len);

Vt1 = zeros(1, t_len);
Vt2 = zeros(1, t_len);
Vt3 = zeros(1, t_len);
Vt4 = zeros(1, t_len);
Vt5 = zeros(1, t_len);
Vt6 = zeros(1, t_len);

da = zeros(1, t_len); % duty cycle array for phase leg a
db = zeros(1, t_len); % duty cycle array for phase leg b
dc = zeros(1, t_len); % duty cycle array for phase leg c

%% while loop for simulation
k = 1; % while loop counter
while (t(k) < tend)
    %% Phases and Angles 
    theta_ac(k + 1) = (w_ac * t(k)) * 180 / pi;
    theta_ac(k + 1) = mod(theta_ac(k + 1), 360);
    %% Duty Cycles
    da(k) = m * cosd(theta_ac(k) - phi);
    db(k) = m * cosd(theta_ac(k) - phi - phi_dis);
    dc(k) = m * cosd(theta_ac(k) - phi + phi_dis);
    %% Switching logic
    if (da(k) >= tri(k))
        T1 = 1;
        T4 = 0;
    elseif (da(k) < tri(k))
        T1 = 0;
        T4 = 1;
    end
    if (db(k) >= tri(k))
        T2 = 1;
        T5 = 0;
    elseif (db(k) < tri(k))
        T2 = 0;
        T5 = 1;
    end
    if (dc(k) >= tri(k))
        T3 = 1;
        T6 = 0;
    elseif (dc(k) < tri(k))
        T3 = 0;
        T6 = 1;
    end
    %% Phase Leg, Diodes, and Transistors Voltages
    if T1 == 1
        Vag(k) = Vdc;
        Vt1(k) = 0;
        Vd1(k) = 0;
        Vt4(k) = Vdc;
        Vd4(k) = -Vdc;
    elseif T4 == 1
        Vag(k) = 0;
        Vt4(k) = 0;
        Vd4(k) = 0;
        Vt1(k) = Vdc;
        Vd1(k) = -Vdc;
    end
    if T2 == 1
        Vbg(k) = Vdc;
        Vt2(k) = 0;
        Vd2(k) = 0;
        Vt5(k) = Vdc;
        Vd5(k) = -Vdc;
    elseif T5 == 1
        Vbg(k) = 0;
        Vt5(k) = 0;
        Vd5(k) = 0;
        Vt2(k) = Vdc;
        Vd2(k) = -Vdc;
    end
    if T3 == 1
        Vcg(k) = Vdc;
        Vt3(k) = 0;
        Vd3(k) = 0;
        Vt6(k) = Vdc;
        Vd6(k) = -Vdc;
    elseif T6 == 1
        Vcg(k) = 0;
        Vt6(k) = 0;
        Vd6(k) = 0;
        Vt3(k) = Vdc;
        Vd3(k) = -Vdc;
    end
    %% AC Voltages and Currents
    Vas(k) = (2 / 3) * Vag(k) - (1 / 3) * Vbg(k) - (1 / 3) * Vcg(k);
    Vbs(k) = (2 / 3) * Vbg(k) - (1 / 3) * Vag(k) - (1 / 3) * Vcg(k);
    Vcs(k) = (2 / 3) * Vcg(k) - (1 / 3) * Vbg(k) - (1 / 3) * Vag(k);
    
    Ias(k + 1) = Ias(k) + dt * (Vas(k) - r * Ias(k)) / l;
    Ibs(k + 1) = Ibs(k) + dt * (Vbs(k) - r * Ibs(k)) / l;
    Ics(k + 1) = Ics(k) + dt * (Vcs(k) - r * Ics(k)) / l;
    %% Diodes and Transistors Currents
    if (T1 == 1)
        if (Ias(k) >= 0) 
            It1(k) = Ias(k);
            Id1(k) = 0;
        else
            Id1(k) = -Ias(k);
            It1(k) = 0;
        end
    end
    if (T2 == 1)
        if (Ibs(k) >= 0) 
            It2(k) = Ibs(k);
            Id2(k) = 0;
        else
            Id2(k) = -Ibs(k);
            It2(k) = 0;
        end
    end
    if (T3 == 1)
        if (Ics(k) >= 0) 
            It3(k) = Ics(k);
            Id3(k) = 0;
        else
            Id3(k) = -Ics(k);
            It3(k) = 0;
        end
    end
    if (T4 == 1)
        if (Ias(k) >= 0) 
            Id4(k) = Ias(k);
            It4(k) = 0;
        else
            It4(k) = -Ias(k);
            Id4(k) = 0;
        end
    end
    if (T5 == 1)
        if (Ibs(k) >= 0) 
            Id5(k) = Ibs(k);
            It5(k) = 0;
        else
            It5(k) = -Ibs(k);
            Id5(k) = 0;
        end
    end
    if (T6 == 1)
        if (Ics(k) >= 0) 
            Id6(k) = Ics(k);
            It6(k) = 0;
        else
            It6(k) = -Ics(k);
            Id6(k) = 0;
        end
    end

    Idc(k) = It1(k) - Id1(k) + It2(k) - Id2(k) + It3(k) - Id3(k);
    % Time array and steps update
    t(k + 1) = t(k) + dt;
    k = k + 1;
end

%% Calculate amplitudes of fundamental frequency component of Vas and Ias
[avg, ak, bk, rcon, err] = fourseries(t, Vas, T, 100);
Vamp = sqrt((ak(1)) ^ 2 + (bk(1)) ^ 2);
disp(Vamp);
[avg, ak, bk, rcon, err] = fourseries(t, Ias, T, 100);
Iamp = sqrt((ak(1)) ^ 2 + (bk(1)) ^ 2);
disp(Iamp);

%% Plots
figure(1)
subplot(2,1,1)
plot(theta_ac(200003:end), Vas(200003:end), ...
    theta_ac(200003:end), Vbs(200003:end), ...
    theta_ac(200003:end), Vcs(200003:end));
title('AC Voltages(V) vs \theta_{AC}(degree)')
grid on
legend('Vas','Vbs','Vcs')
subplot(2,1,2)
plot(theta_ac(200003:end), Ias(200003:end), ...
    theta_ac(200003:end), Ibs(200003:end), ...
    theta_ac(200003:end), Ics(200003:end));
title('AC Currents(A) vs \theta_{AC}(degree)')
grid on
legend('Ias','Ibs','Ics')

figure(2)
subplot(3,2,1)
plot(theta_ac(200003:end), Vd1(200003:end));
%xlim([0 360])
title('Vd1(V) vs \theta_{AC}(degree)')
subplot(3,2,3)
plot(theta_ac(200003:end), Vd2(200003:end));
%xlim([0 360])
title('Vd2(V) vs \theta_{AC}(degree)')
subplot(3,2,5)
plot(theta_ac(200003:end), Vd3(200003:end));
%xlim([0 360])
title('Vd3(V) vs \theta_{AC}(degree)')
subplot(3,2,2)
plot(theta_ac(200003:end), Vd4(200003:end));
%xlim([0 360])
title('Vd4(V) vs \theta_{AC}(degree)')
subplot(3,2,4)
plot(theta_ac(200003:end), Vd5(200003:end));
%xlim([0 360])
title('Vd5(V) vs \theta_{AC}(degree)')
subplot(3,2,6)
plot(theta_ac(200003:end), Vd6(200003:end));
%xlim([0 360])
title('Vd6(V) vs \theta_{AC}(degree)')

figure(3)
subplot(3,2,1)
plot(theta_ac(200003:end), Vt1(200003:end));
%xlim([0 360])
title('Vt1(V) vs \theta_{AC}(degree)')
subplot(3,2,3)
plot(theta_ac(200003:end), Vt2(200003:end));
%xlim([0 360])
title('Vt2(V) vs \theta_{AC}(degree)')
subplot(3,2,5)
plot(theta_ac(200003:end), Vt3(200003:end));
%xlim([0 360])
title('Vt3(V) vs \theta_{AC}(degree)')
subplot(3,2,2)
plot(theta_ac(200003:end), Vt4(200003:end));
%xlim([0 360])
title('Vt4(V) vs \theta_{AC}(degree)')
subplot(3,2,4)
plot(theta_ac(200003:end), Vt5(200003:end));
%xlim([0 360])
title('Vt5(V) vs \theta_{AC}(degree)')
subplot(3,2,6)
plot(theta_ac(200003:end), Vt6(200003:end));
%xlim([0 360])
title('Vt6(V) vs \theta_{AC}(degree)')

figure(4)
subplot(3,2,1)
plot(theta_ac(200003:end), Id1(200003:end));
%xlim([0 360])
title('Id1(A) vs \theta_{AC}(degree)')
subplot(3,2,3)
plot(theta_ac(200003:end), Id2(200003:end));
%xlim([0 360])
title('Id2(A) vs \theta_{AC}(degree)')
subplot(3,2,5)
plot(theta_ac(200003:end), Id3(200003:end));
%xlim([0 360])
title('Id3(A) vs \theta_{AC}(degree)')
subplot(3,2,2)
plot(theta_ac(200003:end), Id4(200003:end));
%xlim([0 360])
title('Id4(A) vs \theta_{AC}(degree)')
subplot(3,2,4)
plot(theta_ac(200003:end), Id5(200003:end));
%xlim([0 360])
title('Id5(A) vs \theta_{AC}(degree)')
subplot(3,2,6)
plot(theta_ac(200003:end), Id6(200003:end));
%xlim([0 360])
title('Id6(A) vs \theta_{AC}(degree)')

figure(5)
subplot(3,2,1)
plot(theta_ac(200003:end), It1(200003:end));
%xlim([0 360])
title('It1(A) vs \theta_{AC}(degree)')
subplot(3,2,3)
plot(theta_ac(200003:end), It2(200003:end));
%xlim([0 360])
title('It2(A) vs \theta_{AC}(degree)')
subplot(3,2,5)
plot(theta_ac(200003:end), It3(200003:end));
%xlim([0 360])
title('It3(A) vs \theta_{AC}(degree)')
subplot(3,2,2)
plot(theta_ac(200003:end), It4(200003:end));
%xlim([0 360])
title('It4(A) vs \theta_{AC}(degree)')
subplot(3,2,4)
plot(theta_ac(200003:end), It5(200003:end));
%xlim([0 360])
title('It5(A) vs \theta_{AC}(degree)')
subplot(3,2,6)
plot(theta_ac(200003:end), It6(200003:end));
%xlim([0 360])
title('It6(A) vs \theta_{AC}(degree)')

figure(6)
plot(theta_ac(200003:end), Idc(200003:end));
%xlim([0 360])
title('Idc(A) vs \theta_{AC}(degree)')

figure(7)
plot(theta_ac(200003:end), Vag(200003:end), theta_ac(200003:end), ...
    Vbg(200003:end), theta_ac(200003:end), Vcg(200003:end));
title('Phase leg Voltages(V) vs \theta_{AC}(degree)')
grid on
legend('Vag','Vbg','Vcg')

figure(8)
plot(t(1:end-1), da(1:end-1), t(1:end-1), db(1:end-1), t(1:end-1), dc(1:end-1));
title('Duty Cycles vs \theta_{AC}(degree)')
grid on
legend('da','db','dc')
