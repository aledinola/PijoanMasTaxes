function K_to_L = fun_KL_from_r(r, theta, delta)
%FUN_KL_FROM_R Capital-labor ratio implied by the firm's interest-rate FOC.
%
% Inputs
%   r                      : net interest rate
%   theta,delta            : production labor share and depreciation rate
%
% Output
%   K_to_L                 : capital-labor ratio consistent with r

K_to_L = ((r + delta) / (1 - theta))^(-1 / theta);

end %end function
