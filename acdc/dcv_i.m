%% dcv_i.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the average output voltage versus
current for modes 1, 2, and 3. It takes in E_amp - back EMF
amplitude, x_ac - product of fundamental angular frequency and AC source 
inductance. It returns the modes' voltages and currents to plot.dc
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[V_m1, i_m1, V_m2, i_m2, V_m3, i_m3] = dcv_i(E_amp, x_ac)
    %% Mode 1 Average Voltage and DC Current
    for gamma = 0:60
        V_m1(gamma + 1) = 3 * sqrt(3) * E_amp * (1 + cosd(gamma)) / (2 * pi);
        i_m1(gamma + 1) = sqrt(3) * E_amp * (1 - cosd(gamma)) / (2 * x_ac);
    end

    %% Mode 2 Average Voltage and DC Current
    for alpha = 0:30
        V_m2(alpha + 1) = 9 * E_amp * cosd(alpha + 30) / (2 * pi);
        i_m2(alpha + 1) = sqrt(3) * E_amp * sind(alpha + 30) / (2 * x_ac);
    end

    %% Mode 3 Average Voltage and DC Current
    for delta = 0:60
        V_m3(delta + 1) = 9 * E_amp * (1 - sind(delta + 30)) / (2 * pi);
        i_m3(delta + 1) = E_amp * (1 + sind(delta + 30)) / (2 * x_ac);
    end