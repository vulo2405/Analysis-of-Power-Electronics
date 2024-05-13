%% average.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dt is time step, x is the waveform that is being analyzed, T is its fundamental 
period of the waveform.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[avg] = average(dt,x,T)
    %% Parameters
    avg = 0;
    N = T / dt; % Number of subintervals
    k = length(x);
    %% Calculate average (avg) of a waveform using Riemann sum 
    for i = k:-1:k-N
        avg = avg + dt * x(i);
    end
    avg = avg / T;