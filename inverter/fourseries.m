%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
t is time, x is the waveform that is being analyzed, T is its fundamental 
period of the waveform, N is the number of terms of the fourier series 
that is desired.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fourseries.m
function[avg,ak,bk,rcon,err] = fourseries(t,x,T,N)
    %% Parameters
    avg = 0;
    err = 0;
    dt = t(2) - t(1); % Size of subinterval
    sample = T / dt; % Number of subintervals
    ak = zeros(1, N);
    bk = zeros(1, N);
    w = 2 * pi / T;
    rcon = zeros(size(t));
    
    %% Calculate average (avg) of a waveform using Riemann sum 
    for i = 1:sample
        avg = avg + dt * x(i);
    end
    avg = avg / T;
    
    %% Calculate ak and bk using Riemann sum
    for k = 1:N
        for n = 1:sample
            ak(k) = ak(k) + x(n) * dt * cos(k * w * t(n));
            bk(k) = bk(k) + x(n) * dt * sin(k * w * t(n));
        end
    end
    ak = 2 * ak / T;
    bk = 2 * bk / T;
    
    %% Calculate rcon, the reconstructed approximation of x
    % First we get a0 term, which is also the average
    rcon = avg;
    % Reconstruct x
    for k = 1:N
        rcon = rcon + ak(k) * cos(k * w .* t) + bk(k) * sin(k * w .* t);
    end

    %% Calculate error
    for i = 1:sample
        err = err + (rcon(i) -  x(i)) ^ 2;
    end
    err = sqrt(err / sample);
end