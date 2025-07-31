function F = Model_ReturnFn(d_val,aprime,a,z,K_to_L,crra,lambda,nu,theta,delta)
% The return function is essentially the combination of the utility
% function and the constraints.

[r,w] = fun_prices(K_to_L,theta,delta);

income = w*z*d_val+r*a;

c = income+a-aprime; % Budget Constraint

F = -inf;

if c>0
    % WARNING: this will not work if crra=1 and/or nu=1
    F = (c^(1-crra)-1)/(1-crra);
    F = F + lambda*((1-d_val)^(1-nu )-1)/(1-nu);
end % end if

end %end function