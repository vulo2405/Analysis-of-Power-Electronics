%% genproc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the system and plot the required graphs.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%% Initial Parameters
E_amp = 392; 
w = 377;
T = 2 * pi / w;
phi = 120; % Phase angle
l = 1e-3;
x_ac = w * l;
gam1 = 45;
gam2 = 60;
d = 15;
l_dc = 5e-3;
%% Time and Angle Arrays
dt = 1e-7;
t = 0:dt:T;
t2 = 0:dt:4*T; % Same in forwardEuler.m
theta = 0:0.01:360;
theta_len = length(theta);
%% Arrays of Mode Voltages and Currents 
V_m1 = zeros(1,60);
V_m2 = zeros(1,30);
V_m3 = zeros(1,60);
i_m1 = zeros(1,60);
i_m2 = zeros(1,30);
i_m3 = zeros(1,60);
%% Arrays of Phase Currents in Commutation and Delay Cases
i_as45 = zeros(1,theta_len);
i_bs45 = zeros(1,theta_len);
i_cs45 = zeros(1,theta_len);
i_asd = zeros(1,theta_len);
i_bsd = zeros(1,theta_len);
i_csd = zeros(1,theta_len);
%% Arrays of Forward Euler Mode Voltages and Currents
V_m1fe = zeros(1,60);
V_m2fe = zeros(1,30);
V_m3fe = zeros(1,60);
i_m1fe = zeros(1,60);
i_m2fe = zeros(1,30);
i_m3fe = zeros(1,60);
%% Resistance Values for Forward Euler calculations
R_1 = [inf 188.85 47.03 20.77 11.58 7.33 5.01 3.62 2.72 2.1 1.66 1.33 1.08];
R_2 = [1.08 0.89 0.74 0.62 0.52 0.44 0.36];
R = [R_1, R_2];
Rlen = length(R);
Vfe_fe = zeros(1,Rlen); % Voltages calculated using Forward Euler
Ife_fe = zeros(1,Rlen); % Currents calculated using Forward Euler
%% PLOTS & CALCULATIONS
%% Back EMFs vs Theta
e_as = E_amp * cosd(theta);
e_bs = E_amp * cosd(theta - phi);
e_cs = E_amp * cosd(theta + phi);
figure(1)
plot(theta, e_as, theta, e_bs, theta, e_cs)
title('Three Phase Back EMFs vs \theta_{ac}');
legend('Phase a','Phase b','Phase c')
grid on
%% Average Voltages and DC Currents Plot
[V_m1, i_m1, V_m2, i_m2, V_m3, i_m3] = dcv_i(E_amp, x_ac);
figure(2)
plot(i_m1, V_m1, i_m2, V_m2, i_m3, V_m3)
title('Modes Voltages and Currents');
xlabel('Currents (A)')
ylabel('Voltages (V)')
grid on
legend('Mode 1', 'Mode 2', 'Mode3')
%% Phase Currents vs Theta_ac with commutation of 45
[i_as45, i_bs45, i_cs45] = commutation(gam1, x_ac, i_m1(46), theta);
figure(3)
plot(theta, i_as45, theta, i_bs45, theta, i_cs45)
title('Phase Currents with Commutative Load Resistance');
grid on
legend('Phase a', 'Phase b', 'Phase c')
%% Phase Currents vs Theta_ac with commutation of 60 & delay of 15
[i_asd, i_bsd, i_csd] = delay(gam2, x_ac, i_m2(16), theta, d);
figure(4)
plot(theta, i_asd, theta, i_bsd, theta, i_csd)
title('Phase Currents with Commutative and Delay Load Resistance');
grid on
legend('Phase a', 'Phase b', 'Phase c')
%% Maximum Power Load
P = zeros(1,150); % Power array
V = [V_m1, V_m2, V_m3];
I = [i_m1, i_m2, i_m3];
P = V .* I;
P_max = max(P);
P_pos = find(P == P_max);
R_min = V(P_pos) ^ 2 / P_max;
disp(P_max);
disp(R_min);
%% Forward Euler Mode Voltages and Currents
% Note: the indexing seen here in every arrays is to compensate for the mod
% function used earlier to calculate and update theta. This function,
% after hitting 360, returns straight back to zero, which results in an
% overlap of many periods, which affects the integrity of the analysis.
[V_m1fe, i_m1fe, V_m2fe, i_m2fe, V_m3fe, i_m3fe] = dcv_i(E_amp, x_ac);
Vfe_a = [V_m1fe, V_m2fe, V_m3fe];
Ife_a = [i_m1fe, i_m2fe, i_m3fe];
for i = 1:Rlen
    [~, ~, ~, ~, Idc_fe, Vdc_fe] = forwardEuler(R(i), E_amp, T);
    Vdc_fe_fe(i) = average(dt, Vdc_fe, T);
    Idc_fe_fe(i) = average(dt, Idc_fe, T);
end
figure(5)
plot(Idc_fe_fe, Vdc_fe_fe, Ife_a, Vfe_a)
title('Forward Euler and Analysis Mode Voltages and Currents');
xlabel('Currents (A)')
ylabel('Voltages (V)')
grid on
legend('Forward Euler', 'Analysis')
%% Commutation Degree of 45 Scenario with Forward Euler Algorithm
Rcomm = 2.098; 
[Ias45, Ibs45, Ics45, theta45, Idc45, ~] = forwardEuler(Rcomm, E_amp, T);
Idc_comm = average(dt, Idc45, T);
disp(Idc_comm);
figure(6) 
plot(t2, Idc45(1:end-1))
yline(Idc_comm,'k')
text(0, Idc_comm, "<i_{DC}>", "HorizontalAlignment", "right");
title('DC Currents and DC Current Average Line for Commutative Angle');
xlabel('Time (s)')
grid on
figure(7)
plot(theta45(500002:end), Ias45(500002:end), ...
     theta45(500002:end), Ibs45(500002:end), ...
     theta45(500002:end), Ics45(500002:end))
xlabel('\theta_{ac} (degrees)')
title('Phase Currents with Commutative Load Resistance with Forward Euler');
grid on
legend('Phase a', 'Phase b', 'Phase c')
%% Commutation Degree and Delay Scenario with Forward Euler Algorithm
Rdelay = 0.62355;
[Ias60, Ibs60, Ics60, theta60, Idc60, ~] = forwardEuler(Rdelay, E_amp, T);
Idc_delay = average(dt, Idc60, T);
disp(Idc_delay);
figure(8) 
plot(t2, Idc60(1:end-1))
yline(Idc_delay,'k')
text(0, Idc_delay, "<i_{DC}>", "HorizontalAlignment", "right");
title('DC Currents and DC Current Average Line for Commutative and Delay Angle');
xlabel('Time (s)')
grid on
figure(9)
plot(theta60(500002:end), Ias60(500002:end), ...
     theta60(500002:end), Ibs60(500002:end), ...
     theta60(500002:end), Ics60(500002:end))
xlabel('\theta_{ac} (degrees)')
title('Phase Currents with Commutative and Delay Load Resistance with Forward Euler');
grid on
legend('Phase a', 'Phase b', 'Phase c')