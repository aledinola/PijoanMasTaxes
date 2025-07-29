function c = Model_cons(d_val,aprime,a,z,K_to_L,theta,delta,lambda_hsv,tau_hsv)
% The return function is essentially the combination of the utility
% function and the constraints.


r = (1-theta)*K_to_L^(-theta)-delta;
w = theta*K_to_L^(1-theta);

income = w*z.*d_val+r*a;
taxes  = income-lambda_hsv*income.^(1-tau_hsv);

c = income+a-taxes-aprime; % Budget Constraint