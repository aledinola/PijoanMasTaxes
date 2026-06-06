function F = Model_ReturnFn(d_val, aprime, a, z, r, crra, lambda, nu, theta, delta)
%MODEL_RETURNFN Household period utility using prices implied by r.
%
% Inputs
%   d_val,aprime,a,z       : hours, next assets, current assets, and productivity
%   r                      : net interest rate
%   crra,lambda,nu         : preference parameters
%   theta,delta            : production labor share and depreciation rate
%
% Output
%   F                      : period utility, or -Inf if consumption is infeasible

c = Model_cons(d_val, aprime, a, z, r, theta, delta);

F = -inf;

if c > 0
    % WARNING: this will not work if crra=1 and/or nu=1.
    F = (c^(1 - crra) - 1) / (1 - crra) + ...
        lambda * ((1 - d_val)^(1 - nu) - 1) / (1 - nu);
end % end if

end %end function
