function flag=valid(x)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


if iscell(x)
    
    flag=all(cellfun(@(x)utils.error.valid(x),x));
    
elseif isa(x,'double')
    
    flag=isreal(x) && ~any(isnan(x(:))) && ~any(isinf(x(:)));
    
elseif isa(x,'tsparse')
    
    flag=utils.error.valid(x.v);
    
else
    
    error(['class ',class(x),' not a valid type for checking validity'])
    
end

% valid=@(x)~any(isnan(x(:))) && ~any(isinf(x(:))); % nans in jacobian
end