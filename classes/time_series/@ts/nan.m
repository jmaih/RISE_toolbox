function db=nan(start_date,varargin)
% INTERNAL FUNCTION
%

data=nan(varargin{:});

db=ts(start_date,data,[],[],true);

