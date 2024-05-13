%% delay.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the system with a delay and a commutation
angle.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[i_as, i_bs, i_cs] = delay(g, x, i_dc, theta_ac, d)
    for k = 1:length(theta_ac)
        %% Phase a current
        if(theta_ac(k) >= 0 && theta_ac(k) <= d)
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) + 120) - sind(315) ...
                      + sind(435)) / x;
        elseif(theta_ac(k) > d && theta_ac(k) < (g + d))
            i_as(k) = i_dc;
        elseif(theta_ac(k) >= (g + d) && theta_ac(k) < (60 + g + d))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) - 120) - sind(g + d) ...
                      + sind(g + d - 120)) / x + i_dc;
        elseif(theta_ac(k) >= (60 + g + d) && theta_ac(k) < (120 + g + d))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) + 120) - sind(60 + g + d) ...
                      + sind(g + d + 180)) / x;
        elseif(theta_ac(k) >= (120 + g + d) && theta_ac(k) < (180 + g + d))
            i_as(k) = -i_dc;
        elseif(theta_ac(k) >= (180 + g + d) && theta_ac(k) < (240 + g + d))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) - 120) - sind(180 + g +d) ...
                      + sind(g + d + 60)) / x - i_dc;
         elseif(theta_ac(k) >= (240 + g + d) && theta_ac(k) <= 360)
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) + 120) - sind(240 + g + d) ...
                      + sind(g + d + 360)) / x;
        end
    end
    disp(196/x);
    shift_num = round(length(theta_ac) * 120 / 360);
    i_bs = circshift(i_as, [0, shift_num]);
    i_cs = circshift(i_as, [0, -shift_num]);
end
    