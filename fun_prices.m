function [r,w] = fun_prices(K_to_L,theta,delta)
%FUN_PRICES Firm FOC prices for a Cobb-Douglas production function.
%
% Inputs
%   K_to_L                 : capital-labor ratio
%   theta,delta            : production labor share and depreciation rate
%
% Output
%   r,w                    : net interest rate and wage

r = (1-theta)*K_to_L^(-theta)-delta;
w = theta*K_to_L^(1-theta);

end %end function
