function c = Model_cons(d_val,aprime,a,z,K_to_L,theta,delta)
% The return function is essentially the combination of the utility
% function and the constraints.

[r,w] = fun_prices(K_to_L,theta,delta);

income = w*z.*d_val+r*a;

c = income+a-aprime; % Budget Constraint

end