function db=rand(start_date,varargin)
% INTERNAL FUNCTION
%

data=rand(varargin{:});

db=ts(start_date,data);

