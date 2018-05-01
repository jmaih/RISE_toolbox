function this3=and(this1,this2)
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
