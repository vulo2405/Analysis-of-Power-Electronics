%% provefourseries.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This is the script to call the function fourseries to prove its correctness. 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

T = 1 / 40; % Fundamental Period
t = 0:1e-7:5*T; % Time array
 x = 15 - 30*cos(2 * pi * 40 * t) + 20*sin(2 * pi * 40 * t)...
     + 8*cos(2 * pi * 800 * t) + 2*sin(2 * pi * 1600 * t)...
     - 5*cos(2 * pi * 2000 * t); 
N = 50; % Number of Fourier series terms

[avg, ak, bk, rcon, err] = fourseries(t, x, T, N);

disp(avg);
disp(err);

figure(1)
plot(t, x, t, rcon)
xlabel('time(s)')
legend('og wave', 'rcon')
grid on

