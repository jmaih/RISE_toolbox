function this3=mtimes(this1,this2)
% Overloaded times function for ts object
%
% Note:
%    Matrix product does not make sense for time series object, so
%    automatically defaults to element-wise product 
%

this3=times(this1,this2);
end
