%buck2.m
k = 1;
while (t(k) < tend)
    state = sw(D, t(k));
    I_L(k + 1) = (Vin * state - Vload(k) - ...
                 VT - RT * I_L(k)) * (dt / L) + I_L(k);
    I_C(k + 1) = I_L(k) - Vload(k) / Rload;
    Vload(k + 1) = Vload(k) + (dt / C) * (I_L(k) - Vload(k) / Rload);
    I_sw(k) = I_L(k) * state;
    I_D(k) = I_L(k) * (1 - state);
    P_D(k) = state * (1 + I_D(k) * RD) + ...
                (1 - state) * (VD - Vin - I_sw(k) * RD) * (I_D(k));
    P_T(k) = state * (1 + I_sw(k) * RT) + ...
                (1 - state) * (Vin - VT - I_D(k) * RT) * (I_sw(k));
    t(k + 1) = t(k) + dt;
    k = k + 1;
end
