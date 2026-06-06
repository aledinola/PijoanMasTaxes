function corr_xy = fun_corr(x,y,prob)
%FUN_CORR Weighted correlation between two arrays.
%
% Inputs
%   x,y                    : equally sized arrays to correlate
%   prob                   : probability mass over the same entries
%
% Output
%   corr_xy                : weighted correlation between x and y

x = x(:);
y = y(:);
prob = prob(:);

if numel(x)~=numel(y)
    error('Inputs x and y in fun_corr have incompatible size!')
end

% Weighted means
mean_x   = sum(prob.*x);
mean_y   = sum(prob.*y);

% Weighted standard deviation
stddev_x = sqrt(sum(prob.*(x-mean_x).^2));
stddev_y = sqrt(sum(prob.*(y-mean_y).^2));

% Weighted covariance
cov_xy  = sum(prob.*x.*y)-mean_x*mean_y;

% Weighted correlation
corr_xy = cov_xy/(stddev_x*stddev_y);

end %end function
