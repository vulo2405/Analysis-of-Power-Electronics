%buck.m
k = 1;
while (t(k) < tend)
    state = sw(D, t(k));
    I_L(k + 1) = (Vin * state - Vload(k)) * dt / L + I_L(k);
    I_C(k + 1) = I_L(k) - Vload(k) / Rload;
    Vload(k + 1) = Vload(k) + dt * (I_L(k) - Vload(k) / Rload) / C;
    I_sw(k) = I_L(k) * state;
    I_D(k) = I_L(k) * (1 - state);
    t(k + 1) = t(k) + dt;
    k = k + 1;
end
