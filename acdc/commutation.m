%% commutation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
This function is to simulate the system with a commutation angle.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[i_as, i_bs, i_cs] = commutation(g, x, i_dc, theta_ac)
    for k = 1:length(theta_ac)
        %% Phase a current
        if(theta_ac(k) >= 0 && theta_ac(k) <= 60)
            i_as(k) = i_dc;
        elseif(theta_ac(k) > 60 && theta_ac(k) < (60 + g))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) - 120) ...
                      - sqrt(3)) / x + i_dc;
        elseif(theta_ac(k) >= (60 + g) && theta_ac(k) < (120))
            i_as(k) = 0;
        elseif(theta_ac(k) >= 120 && theta_ac(k) < (120 + g))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) + 120) ...
                      - sqrt(3)) / x;
        elseif(theta_ac(k) >= 165 && theta_ac(k) < 240)
            i_as(k) = -i_dc;
        elseif(theta_ac(k) >= 240 && theta_ac(k) < (240 + g))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) - 120) ...
                      + sqrt(3)) / x - i_dc;
        elseif(theta_ac(k) >= (240 + g) && theta_ac(k) < 300)
            i_as(k) = 0;
        elseif(theta_ac(k) >= 300 && theta_ac(k) < (300 + g))
            i_as(k) = 196 * (sind(theta_ac(k)) - ...
                      sind(theta_ac(k) + 120) ...
                      + sqrt(3)) / x;
        elseif(theta_ac(k) >= (300 + g) && theta_ac(k) <= 360)
            i_as(k) = i_dc;
        end
    end
    
    shift_num = round(length(theta_ac) * 120 / 360);
    i_bs = circshift(i_as, [0, shift_num]);
    i_cs = circshift(i_as, [0, -shift_num]);
end
    