%% forwardEuler.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the system and plot the required graphs for 
the analytical process.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ias, Ibs, Ics, theta, Idc, Vdc] = forwardEuler(R, E_amp, T)
    %% Initial Parameters and Arrays
    w = 377;      
    phi = 120;       
    tau = 1e-5;                 % Time constant
    Ls = 1e-3;                  % AC Source Inductance
    L_dc = 5e-3;                % DC Inductance
    %% Timing Array
    tend = 4 * T;
    dt = 1e-7;
    t = 0:dt:tend;
    tlen = length(t)+1;
    %% Initial Arrays
    theta = zeros(1,tlen);
    e_as = zeros(1,tlen);   
    e_bs = zeros(1,tlen);
    e_cs = zeros(1,tlen);
    Ias = zeros(1,tlen);
    Ibs = zeros(1,tlen);
    Ics = zeros(1,tlen);
    VLa = zeros(1,tlen);
    VLb = zeros(1,tlen);
    VLc = zeros(1,tlen);
    Va = zeros(1,tlen);
    Vb = zeros(1,tlen);
    Vc = zeros(1,tlen);
    id1 = zeros(1,tlen);
    id3 = zeros(1,tlen);
    id5 = zeros(1,tlen);
    Itransfer = zeros(1,tlen); % Transfer function for Idc
    Vdc = zeros(1,tlen);
    Idc = zeros(1,tlen);
    %% Forward Euler Algorithm
    epsilon = 0.01;
    k = 1;
    while t(k) < tend
        %% Phase and Angle Update
        theta(k + 1) = w * t(k) * 180 / pi;
        theta(k + 1) = mod(theta(k + 1), 360);
        %% Source Voltages
        e_as(k) = E_amp * cosd(theta(k));
        e_bs(k) = E_amp * cosd(theta(k) - phi);
        e_cs(k) = E_amp * cosd(theta(k) + phi);
        %% Forward Euler Calculation
        if Ias(k) > epsilon
            VLa(k) = Vdc(k);
        elseif Ias(k) < -epsilon
            VLa(k) = 0;
        elseif ((Ias(k) > -epsilon) && (Ias(k) < epsilon))
            VLa(k) = ((2 * epsilon) ^ (-1) * Ias(k) + 0.5) * Vdc(k);
        end
        if Ibs(k) > epsilon
            VLb(k) = Vdc(k);
        elseif Ibs(k) < -epsilon
            VLb(k) = 0;
        elseif ((Ibs(k) > -epsilon) && (Ibs(k) < epsilon))
            VLb(k) = ((2 * epsilon) ^ (-1) * Ibs(k) + 0.5) * Vdc(k);
        end
        if Ics(k) > epsilon
            VLc(k) = Vdc(k);
        elseif Ics(k) < -epsilon
            VLc(k) = 0;
        elseif ((Ics(k) > -epsilon) && (Ics(k) < epsilon))
            VLc(k) = ((2 * epsilon) ^ (-1) * Ics(k) + 0.5) * Vdc(k);
        end
        %% Update Phase Voltage Arrays
        Va(k) = (2 * VLa(k) - VLb(k) - VLc(k)) / 3;
        Vb(k) = (2 * VLb(k) - VLc(k) - VLa(k)) / 3;
        Vc(k) = (2 * VLc(k) - VLa(k) - VLb(k)) / 3;
        %% Update Phase Current Arrays
        Ias(k+1) = dt * (e_as(k) - Va(k)) / Ls + Ias(k);
        Ibs(k+1) = dt * (e_bs(k) - Vb(k)) / Ls + Ibs(k);
        Ics(k+1) = dt * (e_cs(k) - Vc(k)) / Ls + Ics(k);
        %% Diode Current Logic
        if Ias(k+1) > 0
            id1(k+1) = Ias(k+1);
        else
            id1(k+1) = 0;
        end
        if Ibs(k+1) > 0
            id3(k+1) = Ibs(k+1);
        else
            id3(k+1) = 0;
        end
        if Ics(k+1) > 0
            id5(k+1) = Ics(k+1);
        else
            id5(k+1) = 0;
        end
        %% Calculate Vdc and Idc
        Idc(k+1) = id1(k+1) + id3(k+1) + id5(k+1);
        Itransfer(k+1) = dt * (R * Idc(k) - Vdc(k)) / tau + Itransfer(k);
        Vdc(k+1) = L_dc * Idc(k+1) / tau + Itransfer(k+1);
        %% Time Update
        t(k+1) = t(k) + dt;
        k = k + 1;
    end
end

