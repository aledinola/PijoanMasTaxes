function [r,w] = fun_prices(K_to_L,theta,delta)

r = (1-theta)*K_to_L^(-theta)-delta;
w = theta*K_to_L^(1-theta);

end