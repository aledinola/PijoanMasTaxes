function c = Model_cons(d_val, aprime, a, z, r, theta, delta)
%MODEL_CONS Consumption implied by choices and prices from r.
%
% Inputs
%   d_val,aprime,a,z       : hours, next assets, current assets, and productivity
%   r                      : net interest rate
%   theta,delta            : production labor share and depreciation rate
%
% Output
%   c                      : current consumption from the household budget constraint

w = fun_w_from_r(r, theta, delta);

income = w * z * d_val + r * a;
c = income + a - aprime; % Budget constraint

end %end function
