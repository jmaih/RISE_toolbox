function n=numel(this,varargin)
% Returns the number of data, i.e., (length of time)x(number of
%   variables)x(number of panels) 
%
% ::
%
%    n = numel(db);
%
% Args:
%
%    db (ts object): time series object
%
% Returns:
%    :
%
%       - n (integer): total number of data points, i.e., (length of
%         time)x(number of variables)x(number of panels) 
%

if isempty(varargin)
    
    n=builtin('numel',this);
    
else
    
    n=1;
    
end

end