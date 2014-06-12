function db=randn(start_date,varargin)

data=randn(varargin{:});

db=ts(start_date,data);

