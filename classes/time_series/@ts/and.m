function this3=and(this1,this2)
% INTERNAL FUNCTION: Combines two databases into a combined database
%
% Note:
%    - It is assumed that the frequency and length of time series are the same between the two databases.
%

% this1= ts(1990,rand(10,1))
% this2= ts('1990Q1',rand(10,1))
% this3= this1 & this2

if nargin~=2
    error([mfilename,':: number of arguments should be 2'])
end
if isempty(this1)
    this3=this2;
elseif isempty(this2)
    this3=this1;
else
    this3=cat(2,this1,this2);
end

end
