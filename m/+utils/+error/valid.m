function flag=valid(x)

flag=~any(isnan(x(:))) && ~any(isinf(x(:)));

% valid=@(x)~any(isnan(x(:))) && ~any(isinf(x(:))); % nans in jacobian
end