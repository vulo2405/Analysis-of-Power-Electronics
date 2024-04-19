%sw.m
function [st] = sw(D, t)
aksum = 1 / 2; % a0 value
for k = 1:30
    aksum = aksum + (2 * cos(pi * k) - cos(2 * pi * k) - 1) / ((pi ^ 2) ...
            * (k ^ 2)) * cos(20000 * pi * k * t);
end
if aksum >= D
    st = 1;
else
    st = 0;
end
end


