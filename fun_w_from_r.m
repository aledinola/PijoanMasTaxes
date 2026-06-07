function [w, K_to_L] = fun_w_from_r(r, theta, delta)
%FUN_W_FROM_R Wage and capital-labor ratio implied by the interest-rate FOC.
%
% Inputs
%   r                      : net interest rate
%   theta,delta            : production labor share and depreciation rate
%
% Outputs
%   w                      : wage implied by r
%   K_to_L                 : capital-labor ratio implied by r

K_to_L = ((r + delta) / (1 - theta))^(-1 / theta);
w = theta * K_to_L^(1 - theta);

end %end function
