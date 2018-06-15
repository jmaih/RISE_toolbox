function db=randn(start_date,varargin)
% INTERNAL FUNCTION
%

data=randn(varargin{:});

db=ts(start_date,data);

