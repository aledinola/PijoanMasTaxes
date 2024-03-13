function [r,w] = fun_prices(KL,par)

theta = par.theta;
delta = par.delta;

r = (1-theta)*KL^(-theta)-delta;
w = theta*KL^(1-theta);

end