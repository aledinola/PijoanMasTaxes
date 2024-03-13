function F = Model_ReturnFn(d_val,aprime,a,z,K_to_L,crra,lambda,nu,theta,delta,lambda_hsv,tau_hsv)
% The return function is essentially the combination of the utility
% function and the constraints.


r = (1-theta)*K_to_L^(-theta)-delta;
w = theta*K_to_L^(1-theta);

income = w*z*d_val+r*a;
taxes  = income-lambda_hsv*income^(1-tau_hsv);

c = income+a-taxes-aprime; % Budget Constraint

F = -inf;

if c>0
    % WARNING: this will not work if crra=1 and/or nu=1
    F = (c^(1-crra)-1)/(1-crra);
    F = F + lambda*((1-d_val)^(1-nu )-1)/(1-nu);
end

end