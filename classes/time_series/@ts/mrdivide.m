function this3=mrdivide(this1,this2)
% Overloaded mrdivide for ts object
%
% Note:
%    It matrix multiplication does not make sense for a time series object, so automatically defaults to element-wise division
%

this3=rdivide(this1,this2);
end
