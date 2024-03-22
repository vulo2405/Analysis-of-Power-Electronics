%buckproc.m
L = 8.8 * 10 ^ -5;
C = 28.409 * 10 ^ -5;
f = 10000;
T = 1 / f;
dt = T / 100;
tend = 80 * T;
t = 0:dt:tend;
t_len = length(t);
t_2 = t_len:-1:(t_len - 2 * T / dt);

D = 0.509;

Vin = 800;
R_heavy = 0.64;

Vload = zeros(1,t_len);
I_L = zeros(1,t_len);
I_C = zeros(1,t_len);
I_sw = zeros(1,t_len);
I_D = zeros(1,t_len);
P_T = zeros(1,t_len);
P_D = zeros(1,t_len);

Vload(1) = 400;
I_L(1) = 506.4;
RT = 0.01;
RD = 0.01;
VT = 1;
VD = 1;

Rload = R_heavy;
buck2;

aVout = aver(Vload, tend, dt);
aIL = aver(I_L, tend, dt);
eff = (aVout ^ 2 / Rload) / (Vin * D * aIL);
disp(aVout);
disp(aIL);
disp(eff);

figure
plot(t_2, P_T(t_2))
title('Transistor Power Loss vs time')
xlabel('time(s)')
ylabel('Transistor Power Loss(W)')

figure
plot(t_2, P_D(t_2))
title('Diode Power Loss vs time')
xlabel('time(s)')
ylabel('Diode Power Loss(W)')

figure
plot(t_2, I_L(t_2))
title('Inductor Current vs time')
xlabel('time(s)')
ylabel('Inductor Current(A)')

figure
plot(t_2, Vload(t_2))
title('Output Voltage vs time')
xlabel('time(s)')
ylabel('Output Voltage(V)')