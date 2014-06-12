function db=rand(start_date,varargin)

data=rand(varargin{:});

db=ts(start_date,data);

