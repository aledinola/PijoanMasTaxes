clear, clc

n = 1000;

% GPU single
x_single = gpuArray.ones(n,1,'single');

% GPU double
x_double = gpuArray.ones(n,1,'double');

% Sum single + double on GPU
x_sum = x_single + x_double;

% Check types
classUnderlying(x_single)
classUnderlying(x_double)
classUnderlying(x_sum)