%aver.m
function [av] = aver(x, T, dt)
sum = 0;
t = 0:dt:T;
t_len = length(t);
dx = T / dt;
for i = t_len:-1:(t_len - dx)
    sum = sum + dt * x(i);
end
av = sum / T;
end