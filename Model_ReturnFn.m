function F = Model_ReturnFn(d_val,aprime,a,z,K_to_L,crra,lambda,nu,theta,delta)

c = Model_cons(d_val,aprime,a,z,K_to_L,theta,delta);

F = -inf;

if c>0
    % WARNING: this will not work if crra=1 and/or nu=1
    F = (c^(1-crra)-1)/(1-crra)+lambda*((1-d_val)^(1-nu )-1)/(1-nu);
end % end if

end %end function